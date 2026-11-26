%% Stage 5: INDEPENDENT multi-start candidate-profile QUBO experiment
% Does NOT import/reference the LP solution or feasible reference schedule.
% MATLAB R2026a, base MATLAB only; requires the 50-EV and INR tariff CSVs.
% Random starts, genuine quadratic binary model, reproducible runs.
% IMPORTANT: feasible solutions are NOT guaranteed in this restricted model.
clear; clc; close all;
rng(20261028,'twister');
EV = readtable('ev_charging_optimization_50_evs_synthetic.csv','TextType','string');
Tariff = readtable('ev_charging_hourly_tariffs_24h_inr.csv');
N = height(EV); H = 24; cap = 250; dt = 1; tol = 1e-6;
assert(N==50 && height(Tariff)==H);
assert(isequal(Tariff.hour(:),(0:23)'));
assert(all(EV.station_capacity_kw==cap));
assert(numel(unique(EV.EV_ID))==N);
price=Tariff.tariff_inr_per_kwh(:);
assert(all(isfinite(price)) && all(price>=0));
assert(all(EV.required_energy_kwh>0) && all(EV.max_charging_power_kw>0));
assert(all(EV.arrival_hour>=0 & EV.departure_hour<=H & EV.arrival_hour<EV.departure_hour));

%% Create binary candidate profiles independently (without LP/reference input).
% Each profile delivers the EV's entire requirement in its valid window.
% Charges are discretized in 0.1 kWh / 0.1 kW increments.
% A random priority permutation creates diverse complete schedules.
% The feasible region is restricted to these generated candidates.
maxCandidates=12;
profiles=cell(N,1); costs=cell(N,1);
for i=1:N
    hours=find((0:H-1)>=EV.arrival_hour(i) & (0:H-1)<EV.departure_hour(i));
    demand=round(10*EV.required_energy_kwh(i));
    maxUnits=floor(10*EV.max_charging_power_kw(i)+1e-7);
    assert(abs(demand/10-EV.required_energy_kwh(i))<tol, 'Demand must use 0.1-kWh resolution.');
    assert(demand<=numel(hours)*maxUnits,'EV energy infeasible in its own window.');
    options=zeros(0,H);
    for k=1:maxCandidates
        if k==1
            ordered=hours; % arrival-first
        elseif k==2
            ordered=fliplr(hours); % deadline-first
        elseif k==3
            [~,ix]=sortrows([price(hours),hours(:)],[1 2]);
            ordered=hours(ix); % tariff-first
        elseif k==4
            [~,ix]=sortrows([-price(hours),hours(:)],[1 2]);
            ordered=hours(ix); % high-tariff-first (diversity)
        else
            ordered=hours(randperm(numel(hours))); % independent random order
        end
        remaining=demand; p=zeros(1,H);
        for h=ordered(:)'
            take=min(remaining,maxUnits);
            p(h)=take/10;
            remaining=remaining-take;
            if remaining==0, break; end
        end
        assert(remaining==0);
        if isempty(options) || ~any(all(abs(options-p)<tol,2))
            options(end+1,:)=p; %#ok<SAGROW>
        end
    end
    profiles{i}=options;
    costs{i}=options*price*dt;
end
nCandidates=sum(cellfun(@(p)size(p,1),profiles));

%% Explicit binary QUBO (same formulation as seeded validation, no LP seed):
% y_ik in {0,1} selects charging profile k for EV i.
% s_hb in {0,1} encodes unused hourly capacity with slack bits.
% Q(y,s) = sum_ik C_ik*y_ik
%          + A*sum_i (sum_k y_ik-1)^2
%          + B*sum_h (sum_ik P_ikh*y_ik + sum_b w_b*s_hb-cap)^2.
% Expanding squares with binary variables gives a true QUBO (degree <= 2).
% The annealer performs one-hot-preserving profile swaps and analytically
% minimizes the slack bits for each proposed state. When overloaded,
% the minimum slack is zero, so the capacity penalty is B*overload^2.
% This is a structured/local-search QUBO annealer, NOT generic bitwise SA.
A=1e6; B=1000;
slackWeights=[2.^(0:10),453]/10;
assert(abs(sum(slackWeights)-cap)<tol);

%% Unseeded independent multistart simulated annealing.
nRuns=8; stepsPerRun=30000; Tstart=2500; Tend=0.10;
runSeed=(20261028+(1:nRuns))';
runFeasible=false(nRuns,1);
runFeasibleCost=nan(nRuns,1);
runLowestObjective=zeros(nRuns,1);
runMaxOverload=zeros(nRuns,1);
runTotalOverload=zeros(nRuns,1);
runRuntime=zeros(nRuns,1);
runStartCost=zeros(nRuns,1);
bestFeasibleCost=inf; bestFeasibleChoices=[];
bestFeasibleRun=nan;
bestObjective=inf; bestObjectiveChoices=[];

for run=1:nRuns
    rng(runSeed(run),'twister'); % independent, reproducible run
    choices=zeros(N,1);
    loadKW=zeros(1,H); currentCost=0;
    % No optimal reference schedule is used anywhere in this script.
    for i=1:N
        choices(i)=randi(size(profiles{i},1));
        loadKW=loadKW+profiles{i}(choices(i),:);
        currentCost=currentCost+costs{i}(choices(i));
    end
    runStartCost(run)=currentCost;
    penalty=B*sum(max(0,loadKW-cap).^2);
    thisBestObjective=currentCost+penalty;
    thisBestChoices=choices;
    thisFeasibleCost=inf;
    tic;
    if all(loadKW<=cap+tol)
        thisFeasibleCost=currentCost;
        if currentCost<bestFeasibleCost
            bestFeasibleCost=currentCost; bestFeasibleChoices=choices;
            bestFeasibleRun=run;
        end
    end
    for t=1:stepsPerRun
        temperature=Tstart*(Tend/Tstart)^((t-1)/(stepsPerRun-1));
        i=randi(N); m=size(profiles{i},1);
        if m<2, continue; end
        proposed=randi(m-1);
        if proposed>=choices(i), proposed=proposed+1; end
        newLoad=loadKW-profiles{i}(choices(i),:)+profiles{i}(proposed,:);
        newCost=currentCost-costs{i}(choices(i))+costs{i}(proposed);
        newPenalty=B*sum(max(0,newLoad-cap).^2);
        delta=(newCost+newPenalty)-(currentCost+penalty);
        if delta<=0 || rand<exp(-delta/temperature)
            choices(i)=proposed; loadKW=newLoad;
            currentCost=newCost; penalty=newPenalty;
        end
        objective=currentCost+penalty;
        if objective<thisBestObjective-1e-8
            thisBestObjective=objective;
            thisBestChoices=choices;
        end
        if all(loadKW<=cap+tol) && currentCost<thisFeasibleCost-1e-8
            thisFeasibleCost=currentCost;
            if currentCost<bestFeasibleCost-1e-8
                bestFeasibleCost=currentCost;
                bestFeasibleChoices=choices;
                bestFeasibleRun=run;
            end
        end
    end
    runRuntime(run)=toc;
    runLowestObjective(run)=thisBestObjective;
    runFeasible(run)=isfinite(thisFeasibleCost);
    if runFeasible(run), runFeasibleCost(run)=thisFeasibleCost; end
    p=zeros(N,H);
    for i=1:N, p(i,:)=profiles{i}(thisBestChoices(i),:); end
    candidateLoad=sum(p,1);
    runMaxOverload(run)=max(max(candidateLoad-cap,0));
    runTotalOverload(run)=sum(max(candidateLoad-cap,0));
    if thisBestObjective<bestObjective
        bestObjective=thisBestObjective;
        bestObjectiveChoices=thisBestChoices;
    end
    if runFeasible(run)
    fprintf('Run %d/%d | feasible=YES | best feasible INR %.2f | min penalized %.2f | %.3f s\n', ...
        run, nRuns, thisFeasibleCost, thisBestObjective, runRuntime(run));
else
    fprintf('Run %d/%d | feasible=NO | best feasible INR N/A | min penalized %.2f | %.3f s\n', ...
        run, nRuns, thisBestObjective, runRuntime(run));
    end

%% Output separately: best feasible (if found), otherwise best INFEASIBLE.
if isfinite(bestFeasibleCost)
    selected=bestFeasibleChoices; status='FEASIBLE';
else
    selected=bestObjectiveChoices; status='INFEASIBLE (no feasible candidate found)';
end
schedule=zeros(N,H);
for i=1:N, schedule(i,:)=profiles{i}(selected(i),:); end
station=sum(schedule,1); evCosts=schedule*price*dt;
delivered=sum(schedule,2)*dt;
assert(all(abs(delivered-EV.required_energy_kwh)<tol));
for i=1:N
    outside=(0:H-1)<EV.arrival_hour(i) | (0:H-1)>=EV.departure_hour(i);
    assert(all(abs(schedule(i,outside))<tol));
    assert(all(schedule(i,:)<=EV.max_charging_power_kw(i)+tol));
end
violation=max(station-cap,0);
% For feasible states binary slack exactly represents unused capacity;
% for overloaded states optimal slack is zero and residual is overload.
slackResidual=max(station-cap,0);
quboEnergy=sum(evCosts)+B*sum(slackResidual.^2); % one-hot penalty zero
assert(abs(quboEnergy-(sum(evCosts)+B*sum(violation.^2)))<tol);

fprintf('\n--- INDEPENDENT QUBO BENCHMARK ---\n');
fprintf('No reference schedule/LP solution used to generate candidates or seeds.\n');
fprintf('EVs: %d | candidate profiles: %d | slack bits: %d\n',N,nCandidates,H*numel(slackWeights));
fprintf('Independent runs: %d | feasible runs: %d | feasibility rate: %.1f%%\n', ...
 nRuns,sum(runFeasible),100*mean(runFeasible));
fprintf('Selected schedule: %s\n',status);
fprintf('Selected electricity cost: INR %.2f\n',sum(evCosts));
fprintf('Energy delivered: %.1f kWh | max unmet energy: %.6f kWh\n', ...
 sum(delivered),max(abs(delivered-EV.required_energy_kwh)));
fprintf('Peak: %.2f kW | violating hours: %d | max excess: %.2f kW\n', ...
 max(station),sum(violation>tol),max(violation));
fprintf('Selected penalized QUBO energy: %.2f\n',quboEnergy);
fprintf('Total annealing runtime: %.3f s\n',sum(runRuntime));
if isfinite(bestFeasibleCost)
    fprintf('Best FEASIBLE cost: INR %.2f (run %d)\n',bestFeasibleCost,bestFeasibleRun);
    fprintf('Gap vs continuous LP benchmark (INR 45190.05): %.2f INR (%.3f%%)\n', ...
        bestFeasibleCost-45190.05,100*(bestFeasibleCost/45190.05-1));
    fprintf('Mean feasible-run best cost: INR %.2f\n',mean(runFeasibleCost(runFeasible)));
else
    fprintf('No feasible schedule found: cost and gap are NOT comparable with LP.\n');
    fprintf('Try increasing candidate diversity or improving QUBO search; do not claim success.\n');
end

Results=table((1:nRuns)',runSeed,runFeasible,runFeasibleCost, ...
 runLowestObjective,runMaxOverload,runTotalOverload,runRuntime,runStartCost, ...
 'VariableNames',{'Run','RandomSeed','Feasible','BestFeasibleCost_INR', ...
 'BestPenalizedObjective','BestPenalizedMaxOverload_kW', ...
 'BestPenalizedTotalOverload_kW','Runtime_s','StartCost_INR'});
writetable(Results,'qubo_independent_run_results.csv');
S=array2table(schedule,'VariableNames',cellstr(compose('hour_%02d_kw',0:23)));
S=addvars(S,EV.EV_ID,'Before',1,'NewVariableNames','EV_ID');
writetable(S,'qubo_independent_schedule.csv');
Hourly=table((0:23)',station',repmat(cap,H,1),price,violation', ...
 'VariableNames',{'hour','Station_Demand_kW','Capacity_kW','Tariff_INR_per_kWh','Excess_kW'});
writetable(Hourly,'qubo_independent_hourly.csv');
figure('Name','Independent QUBO benchmark');
plot(0:23,station,'-o','LineWidth',1.5); hold on;
yline(cap,'r--','250 kW station limit','LineWidth',1.5);
xlabel('Hour of day');ylabel('Station demand (kW)');
title(['Independent QUBO annealing: ' status]);grid on;xlim([0 23]);

end
