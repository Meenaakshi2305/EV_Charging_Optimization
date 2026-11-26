%% Stage 6: Controlled independent QUBO annealing-budget experiment
% MATLAB R2026a; only requires 50-EV input CSV + INR tariff CSV.
% Independent of LP solution and reference schedule; does NOT read either.
% Candidate set, cost/capacity penalty, temperature endpoints held fixed.
% 20 independent starts per budget; same starting choices for each budget.
% Baseline LP value below is for reporting only; NEVER influences search.
% Searches the RESTRICTED candidate-profile QUBO, not all continuous schedules.
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


%% Experiment settings: change only the annealing step budget.
stepBudgets = [5000 15000 30000 60000];
nRuns=20;
runSeeds=(32026000+(1:nRuns))';
Tstart=2500; Tend=0.10;
lpBenchmark=45190.05;  % validated LP cost, reporting only
nB=numel(stepBudgets);

% Generate identical independent random initial profile choices for every
% budget; this allows paired comparisons without LP/reference seeding.
initialChoices=zeros(N,nRuns);
for r=1:nRuns
    rng(runSeeds(r),'twister');
    for i=1:N
        initialChoices(i,r)=randi(size(profiles{i},1));
    end
end

numTests=nRuns*nB;
BudgetSteps=zeros(numTests,1); Run=zeros(numTests,1); RandomSeed=zeros(numTests,1);
Feasible=false(numTests,1); FeasibleCost=nan(numTests,1);
GapPercent=nan(numTests,1); BestPenalized=nan(numTests,1);
BestPenalizedOverload=nan(numTests,1); Runtime_s=nan(numTests,1);
StartCost=nan(numTests,1); StartMaxOverload=nan(numTests,1);

bestOverallCost=inf; bestOverallChoices=[]; bestOverallSteps=nan; bestOverallRun=nan;
experimentTic=tic;
row=0;
for b=1:nB
    steps=stepBudgets(b);
    budgetFeasibleCosts=nan(nRuns,1);
    for r=1:nRuns
        row=row+1;
        BudgetSteps(row)=steps; Run(row)=r; RandomSeed(row)=runSeeds(r);
        rng(runSeeds(r)+500000,'twister');  % annealing stream; independent from candidate generation
        choices=initialChoices(:,r);
        loadKW=zeros(1,H); currentCost=0;
        for i=1:N
            loadKW=loadKW+profiles{i}(choices(i),:);
            currentCost=currentCost+costs{i}(choices(i));
        end
        StartCost(row)=currentCost;
        StartMaxOverload(row)=max(max(loadKW-cap,0));
        penalty=B*sum(max(0,loadKW-cap).^2);
        lowestObjective=currentCost+penalty;
        lowestObjectiveMaxOverload=max(max(loadKW-cap,0));
        runBestCost=inf;
        if all(loadKW<=cap+tol)
            runBestCost=currentCost;
            if currentCost<bestOverallCost-1e-7
                bestOverallCost=currentCost;
                bestOverallChoices=choices;
                bestOverallSteps=steps; bestOverallRun=r;
            end
        end
        runTic=tic;
        for t=1:steps
            temperature=Tstart*(Tend/Tstart)^((t-1)/max(steps-1,1));
            i=randi(N); m=size(profiles{i},1);
            if m<2, continue; end
            proposed=randi(m-1);
            if proposed>=choices(i), proposed=proposed+1; end
            newLoad=loadKW-profiles{i}(choices(i),:)+profiles{i}(proposed,:);
            newCost=currentCost-costs{i}(choices(i))+costs{i}(proposed);
            newPenalty=B*sum(max(0,newLoad-cap).^2);
            delta=(newCost+newPenalty)-(currentCost+penalty);
            if delta<=0 || rand<exp(-delta/temperature)
                choices(i)=proposed;
                loadKW=newLoad;
                currentCost=newCost;
                penalty=newPenalty;
            end
            objective=currentCost+penalty;
            if objective<lowestObjective-1e-8
                lowestObjective=objective;
                lowestObjectiveMaxOverload=max(max(loadKW-cap,0));
            end
            if all(loadKW<=cap+tol) && currentCost<runBestCost-1e-8
                runBestCost=currentCost;
                if currentCost<bestOverallCost-1e-8
                    bestOverallCost=currentCost;
                    bestOverallChoices=choices;
                    bestOverallSteps=steps;
                    bestOverallRun=r;
                end
            end
        end
        Runtime_s(row)=toc(runTic);
        BestPenalized(row)=lowestObjective;
        BestPenalizedOverload(row)=lowestObjectiveMaxOverload;
        Feasible(row)=isfinite(runBestCost);
        if Feasible(row)
            FeasibleCost(row)=runBestCost;
            GapPercent(row)=100*(runBestCost/lpBenchmark-1);
            budgetFeasibleCosts(r)=runBestCost;
        end
    end
    valid=budgetFeasibleCosts(isfinite(budgetFeasibleCosts));
    if isempty(valid)
        fprintf('Budget %6d steps: %2d/%d feasible | no feasible cost | total SA time %.3f s\n', ...
            steps,0,nRuns,sum(Runtime_s(BudgetSteps==steps)));
    else
        fprintf('Budget %6d steps: %2d/%d feasible | best INR %.2f | mean INR %.2f | total SA time %.3f s\n', ...
            steps,numel(valid),nRuns,min(valid),mean(valid),sum(Runtime_s(BudgetSteps==steps)));
    end
end
elapsedTotal=toc(experimentTic);

%% Save all runs, including infeasible ones (NaN costs are intentional).
RunResults=table(BudgetSteps,Run,RandomSeed,Feasible,FeasibleCost,GapPercent, ...
    BestPenalized,BestPenalizedOverload,Runtime_s,StartCost,StartMaxOverload);
writetable(RunResults,'qubo_budget_all_runs.csv');

Budget=zeros(nB,1); FeasibleRuns=zeros(nB,1); FeasibilityRatePct=zeros(nB,1);
BestFeasibleCost_INR=nan(nB,1); MeanFeasibleCost_INR=nan(nB,1);
StdFeasibleCost_INR=nan(nB,1); BestGapPct=nan(nB,1);
MeanGapPct=nan(nB,1); MeanRuntime_s=zeros(nB,1); TotalRuntime_s=zeros(nB,1);
for b=1:nB
    ids=BudgetSteps==stepBudgets(b);
    Budget(b)=stepBudgets(b);
    feasibleCosts=FeasibleCost(ids & Feasible);
    FeasibleRuns(b)=numel(feasibleCosts);
    FeasibilityRatePct(b)=100*FeasibleRuns(b)/nRuns;
    MeanRuntime_s(b)=mean(Runtime_s(ids));
    TotalRuntime_s(b)=sum(Runtime_s(ids));
    if ~isempty(feasibleCosts)
        BestFeasibleCost_INR(b)=min(feasibleCosts);
        MeanFeasibleCost_INR(b)=mean(feasibleCosts);
        StdFeasibleCost_INR(b)=std(feasibleCosts);
        BestGapPct(b)=100*(min(feasibleCosts)/lpBenchmark-1);
        MeanGapPct(b)=100*(mean(feasibleCosts)/lpBenchmark-1);
    end
end
Summary=table(Budget,FeasibleRuns,FeasibilityRatePct,BestFeasibleCost_INR, ...
    MeanFeasibleCost_INR,StdFeasibleCost_INR,BestGapPct,MeanGapPct, ...
    MeanRuntime_s,TotalRuntime_s);
writetable(Summary,'qubo_budget_summary.csv');

fprintf('\n--- CONTROLLED QUBO PERFORMANCE BENCHMARK ---\n');
disp(Summary);
fprintf('SA-only runtime across %d tests: %.3f seconds\n',numTests,sum(Runtime_s));
fprintf('Wall time for full experiment after data/candidate construction: %.3f seconds\n',elapsedTotal);
if isfinite(bestOverallCost)
    fprintf('Best independent feasible cost: INR %.2f at %d steps, run %d\n', ...
        bestOverallCost,bestOverallSteps,bestOverallRun);
    fprintf('LP benchmark: INR %.2f | gap: %.3f%%\n', ...
        lpBenchmark,100*(bestOverallCost/lpBenchmark-1));
    schedule=zeros(N,H);
    for i=1:N
        schedule(i,:)=profiles{i}(bestOverallChoices(i),:);
    end
    station=sum(schedule,1);
    delivered=sum(schedule,2)*dt;
    assert(all(abs(delivered-EV.required_energy_kwh)<tol),'EV energy check failed');
    assert(all(station<=cap+tol),'Capacity check failed');
    for i=1:N
        out=(0:H-1)<EV.arrival_hour(i) | (0:H-1)>=EV.departure_hour(i);
        assert(all(abs(schedule(i,out))<tol),'Charging-window check failed');
        assert(all(schedule(i,:)<=EV.max_charging_power_kw(i)+tol),'EV power check failed');
    end
    assert(abs(sum(sum(schedule,1).*price')-bestOverallCost)<1e-5,'Cost validation failed');
    S=array2table(schedule,'VariableNames',cellstr(compose('hour_%02d_kw',0:23)));
    S=addvars(S,EV.EV_ID,'Before',1,'NewVariableNames','EV_ID');
    writetable(S,'qubo_budget_best_schedule.csv');
    Hourly=table((0:23)',station',repmat(cap,H,1),price, ...
        'VariableNames',{'hour','Station_Demand_kW','Capacity_kW','Tariff_INR_per_kWh'});
    writetable(Hourly,'qubo_budget_best_hourly.csv');
else
    fprintf('No feasible schedule found in any test. Cost/gap cannot be compared.\n');
end

%% Plot cost-quality and runtime trade-off; NaN means no feasible result.
figure('Name','QUBO cost and step budget');
plot(Budget,BestFeasibleCost_INR,'-o','LineWidth',1.5);hold on;
plot(Budget,MeanFeasibleCost_INR,'-s','LineWidth',1.5);
yline(lpBenchmark,'--','LP optimum');grid on;
xlabel('Annealing steps per run');ylabel('Electricity cost (INR)');
title('Independent QUBO: best/mean feasible cost');
legend('Best feasible','Mean of feasible runs','Continuous LP benchmark','Location','best');
figure('Name','QUBO feasibility and runtime');
plot(Budget,FeasibilityRatePct,'-o','LineWidth',1.5);grid on;
xlabel('Annealing steps per run');ylabel('Feasibility rate (%)');
title('Independent QUBO: feasibility rate by search budget');

fprintf('\nSaved: qubo_budget_all_runs.csv, qubo_budget_summary.csv\n');
if isfinite(bestOverallCost)
    fprintf('Saved: qubo_budget_best_schedule.csv, qubo_budget_best_hourly.csv\n');
end
fprintf('Note: only steps vary; all candidate profiles and penalties remain fixed.\n');
fprintf('Mean feasible cost excludes infeasible runs; always report feasibility rate too.\n');
