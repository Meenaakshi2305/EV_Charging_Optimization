%% Stage 7: Independent QUBO candidate-profile sensitivity experiment
% MATLAB R2026a. Same EVs, tariff, capacity and QUBO model as Stage 6.
% No reference schedule or LP output is loaded or used as a starting state.
% IMPORTANT: Only the candidate-profile limit varies. Candidate sets NEST:
% 12 attempted profiles subset of 24 subset of 48. Actual unique counts vary.
% Fixed 30000 SA steps, 20 independent runs, identical penalty / cooling.
clear; clc; close all;
rng(20261028,'twister');
EV=readtable('ev_charging_optimization_50_evs_synthetic.csv','TextType','string');
Tariff=readtable('ev_charging_hourly_tariffs_24h_inr.csv');
N=height(EV); H=24; cap=250; dt=1; tol=1e-6;
assert(N==50 && height(Tariff)==24,'Expected 50 EVs and 24 tariffs.');
assert(isequal(Tariff.hour(:),(0:23)'),'Tariff hours must equal 0:23.');
assert(numel(unique(EV.EV_ID))==N && all(EV.station_capacity_kw==cap));
assert(all(EV.required_energy_kwh>0 & EV.max_charging_power_kw>0));
assert(all(EV.arrival_hour>=0 & EV.departure_hour<=H & EV.arrival_hour<EV.departure_hour));
price=Tariff.tariff_inr_per_kwh(:);
assert(all(isfinite(price) & price>=0));

%% A nested master collection of candidate profiles. No LP/reference data.
% Each candidate supplies ALL energy for one EV in its allowed hours.
% The first 4 attempts use structured hour orderings; others random.
profileLimits=[12 24 48]; maxAttempts=max(profileLimits);
allProfiles=cell(N,1); allCosts=cell(N,1);
% Use separate, fixed stream for profile generation.
rng(20261028,'twister');
for i=1:N
    hours=find((0:H-1)>=EV.arrival_hour(i) & (0:H-1)<EV.departure_hour(i));
    demand=round(10*EV.required_energy_kwh(i));
    maxUnits=floor(10*EV.max_charging_power_kw(i)+1e-7);
    assert(abs(demand/10-EV.required_energy_kwh(i))<tol);
    assert(demand<=numel(hours)*maxUnits,'Individual EV infeasible.');
    attempts=zeros(maxAttempts,H);
    for k=1:maxAttempts
        if k==1
            ordered=hours;
        elseif k==2
            ordered=fliplr(hours);
        elseif k==3
            [~,idx]=sortrows([price(hours),hours(:)],[1 2]);
            ordered=hours(idx);
        elseif k==4
            [~,idx]=sortrows([-price(hours),hours(:)],[1 2]);
            ordered=hours(idx);
        else
            ordered=hours(randperm(numel(hours)));
        end
        remaining=demand; p=zeros(1,H);
        for h=ordered(:)'
            amount=min(maxUnits,remaining);
            p(h)=amount/10; remaining=remaining-amount;
            if remaining==0, break; end
        end
        assert(remaining==0,'Candidate could not deliver EV demand.');
        attempts(k,:)=p;
    end
    allProfiles{i}=attempts;
    allCosts{i}=attempts*price*dt;
end

%% Identical SA settings for each candidate-limit group.
nRuns=20; steps=30000; A=1e6; B=1000;
Tstart=2500; Tend=0.10;
lpBenchmark=45190.05; % ONLY for reporting after search, never used to guide it
nGroups=numel(profileLimits); nTests=nGroups*nRuns;
GroupLimit=zeros(nTests,1); ActualTotalProfiles=zeros(nTests,1);
Run=zeros(nTests,1); RandomSeed=zeros(nTests,1);
Feasible=false(nTests,1); FeasibleCost_INR=nan(nTests,1);
GapPct=nan(nTests,1); BestPenalized_QUBO=nan(nTests,1);
BestPenalizedOverload_kW=nan(nTests,1); Runtime_s=nan(nTests,1);
StartMaxOverload_kW=nan(nTests,1);
runSeeds=(92026000+(1:nRuns))';
% Common random fractions for initialization; indices mapped to each set.
rng(20261029,'twister');
initialUniform=rand(N,nRuns);
bestCost=inf; bestChoice=[]; bestProfiles={}; bestGroup=nan; bestRun=nan;
wallClock=tic; row=0;
for g=1:nGroups
    limit=profileLimits(g);
    profiles=cell(N,1); costs=cell(N,1);
    for i=1:N
        % Remove duplicates while retaining first occurrence / nestedness.
        subset=allProfiles{i}(1:limit,:);
        profiles{i}=unique(subset,'rows','stable');
        costs{i}=profiles{i}*price*dt;
    end
    nProf=sum(cellfun(@(p)size(p,1),profiles));
    fprintf('\nCandidate attempts/EV = %d | unique candidate profiles = %d\n',limit,nProf);
    for r=1:nRuns
        row=row+1;
        GroupLimit(row)=limit; ActualTotalProfiles(row)=nProf;
        Run(row)=r; RandomSeed(row)=runSeeds(r);
        rng(runSeeds(r),'twister');
        choices=zeros(N,1); loadKW=zeros(1,H); currentCost=0;
        for i=1:N
            m=size(profiles{i},1);
            choices(i)=min(m,1+floor(initialUniform(i,r)*m));
            loadKW=loadKW+profiles{i}(choices(i),:);
            currentCost=currentCost+costs{i}(choices(i));
        end
        StartMaxOverload_kW(row)=max(max(loadKW-cap,0));
        penalty=B*sum(max(0,loadKW-cap).^2);
        lowestObjective=currentCost+penalty;
        minPenalizedOverload=max(max(loadKW-cap,0));
        runBestCost=inf;
        if all(loadKW<=cap+tol)
            runBestCost=currentCost;
            if currentCost<bestCost-1e-8
                bestCost=currentCost; bestChoice=choices; bestProfiles=profiles;
                bestGroup=limit; bestRun=r;
            end
        end
        timer=tic;
        for t=1:steps
            temperature=Tstart*(Tend/Tstart)^((t-1)/(steps-1));
            i=randi(N); m=size(profiles{i},1);
            if m<2, continue; end
            proposed=randi(m-1);
            if proposed>=choices(i), proposed=proposed+1; end
            trialLoad=loadKW-profiles{i}(choices(i),:)+profiles{i}(proposed,:);
            trialCost=currentCost-costs{i}(choices(i))+costs{i}(proposed);
            trialPenalty=B*sum(max(0,trialLoad-cap).^2);
            delta=(trialCost+trialPenalty)-(currentCost+penalty);
            if delta<=0 || rand<exp(-delta/temperature)
                choices(i)=proposed; loadKW=trialLoad;
                currentCost=trialCost; penalty=trialPenalty;
            end
            obj=currentCost+penalty;
            if obj<lowestObjective-1e-8
                lowestObjective=obj;
                minPenalizedOverload=max(max(loadKW-cap,0));
            end
            if all(loadKW<=cap+tol) && currentCost<runBestCost-1e-8
                runBestCost=currentCost;
                if currentCost<bestCost-1e-8
                    bestCost=currentCost; bestChoice=choices; bestProfiles=profiles;
                    bestGroup=limit; bestRun=r;
                end
            end
        end
        Runtime_s(row)=toc(timer);
        Feasible(row)=isfinite(runBestCost);
        BestPenalized_QUBO(row)=lowestObjective;
        BestPenalizedOverload_kW(row)=minPenalizedOverload;
        if Feasible(row)
            FeasibleCost_INR(row)=runBestCost;
            GapPct(row)=100*(runBestCost/lpBenchmark-1);
        end
    end
    ids=GroupLimit==limit;
    v=FeasibleCost_INR(ids & Feasible);
    if isempty(v)
        fprintf('Group %2d: %2d/%d feasible | best N/A | mean N/A | SA %.3f s\n', ...
            limit,0,nRuns,sum(Runtime_s(ids)));
    else
        fprintf('Group %2d: %2d/%d feasible | best INR %.2f | mean INR %.2f | SA %.3f s\n', ...
            limit,numel(v),nRuns,min(v),mean(v),sum(Runtime_s(ids)));
    end
end
elapsed=toc(wallClock);

%% Save all 60 runs and their group statistics.
AllRuns=table(GroupLimit,ActualTotalProfiles,Run,RandomSeed,Feasible, ...
    FeasibleCost_INR,GapPct,BestPenalized_QUBO,BestPenalizedOverload_kW, ...
    StartMaxOverload_kW,Runtime_s);
writetable(AllRuns,'qubo_profile_all_runs.csv');
CandidateLimit=zeros(nGroups,1); UniqueProfiles=zeros(nGroups,1);
FeasibleRuns=zeros(nGroups,1); FeasibilityRatePct=zeros(nGroups,1);
BestFeasibleCost_INR=nan(nGroups,1); MeanFeasibleCost_INR=nan(nGroups,1);
StdFeasibleCost_INR=nan(nGroups,1); BestGapPct=nan(nGroups,1);
MeanRuntime_s=zeros(nGroups,1);
for g=1:nGroups
    ids=GroupLimit==profileLimits(g);
    CandidateLimit(g)=profileLimits(g);
    UniqueProfiles(g)=ActualTotalProfiles(find(ids,1));
    good=FeasibleCost_INR(ids & Feasible);
    FeasibleRuns(g)=numel(good);
    FeasibilityRatePct(g)=100*numel(good)/nRuns;
    MeanRuntime_s(g)=mean(Runtime_s(ids));
    if ~isempty(good)
        BestFeasibleCost_INR(g)=min(good);
        MeanFeasibleCost_INR(g)=mean(good);
        StdFeasibleCost_INR(g)=std(good);
        BestGapPct(g)=100*(min(good)/lpBenchmark-1);
    end
end
Summary=table(CandidateLimit,UniqueProfiles,FeasibleRuns,FeasibilityRatePct, ...
    BestFeasibleCost_INR,MeanFeasibleCost_INR,StdFeasibleCost_INR,BestGapPct,MeanRuntime_s);
writetable(Summary,'qubo_profile_summary.csv');
fprintf('\n--- QUBO CANDIDATE-PROFILE SENSITIVITY ---\n');
disp(Summary);
fprintf('SA-only runtime: %.3f s for %d tests\n',sum(Runtime_s),nTests);
fprintf('Wall time (excluding data and initial profile generation): %.3f s\n',elapsed);
if isfinite(bestCost)
    schedule=zeros(N,H);
    for i=1:N
        schedule(i,:)=bestProfiles{i}(bestChoice(i),:);
    end
    demand=sum(schedule,1); energy=sum(schedule,2);
    assert(max(abs(energy-EV.required_energy_kwh))<tol,'Energy validation failed.');
    assert(all(demand<=cap+tol),'Capacity validation failed.');
    for i=1:N
        outside=(0:H-1)<EV.arrival_hour(i) | (0:H-1)>=EV.departure_hour(i);
        assert(all(abs(schedule(i,outside))<tol),'Charging-window check failed.');
        assert(all(schedule(i,:)<=EV.max_charging_power_kw(i)+tol),'EV maximum power check failed.');
    end
    recomputedCost=sum(schedule*price);
    assert(abs(recomputedCost-bestCost)<1e-5,'Cost validation failed.');
    fprintf('Best feasible schedule: INR %.2f | candidate limit %d | run %d\n',bestCost,bestGroup,bestRun);
    fprintf('LP comparison only: INR %.2f | gap %.3f%%\n',lpBenchmark,100*(bestCost/lpBenchmark-1));
    fprintf('Maximum hourly demand %.2f kW | unmet energy %.6f kWh\n',max(demand),max(abs(energy-EV.required_energy_kwh)));
    Out=array2table(schedule,'VariableNames',cellstr(compose('hour_%02d_kw',0:23)));
    Out=addvars(Out,EV.EV_ID,'Before',1,'NewVariableNames','EV_ID');
    writetable(Out,'qubo_profile_best_schedule.csv');
    Hourly=table((0:23)',demand',repmat(cap,H,1),price, ...
        'VariableNames',{'hour','Station_Demand_kW','Capacity_kW','Tariff_INR_per_kWh'});
    writetable(Hourly,'qubo_profile_best_hourly.csv');
else
    fprintf('No feasible schedule: no feasible cost or optimality gap available.\n');
end
figure('Name','QUBO candidate profiles versus cost');
plot(UniqueProfiles,BestFeasibleCost_INR,'-o','LineWidth',1.5); hold on;
plot(UniqueProfiles,MeanFeasibleCost_INR,'-s','LineWidth',1.5);
yline(lpBenchmark,'--','LP benchmark'); grid on;
xlabel('Total unique candidate profiles'); ylabel('Electricity cost (INR)');
title('QUBO profile sensitivity: best and mean feasible costs');
legend('Best feasible','Mean of feasible runs','LP reference','Location','best');
figure('Name','QUBO profile feasibility');
plot(UniqueProfiles,FeasibilityRatePct,'-o','LineWidth',1.5); grid on;
xlabel('Total unique candidate profiles'); ylabel('Feasibility rate (%)');
title('QUBO profile sensitivity: feasibility rate'); ylim([0 105]);
fprintf('Saved: qubo_profile_all_runs.csv, qubo_profile_summary.csv\n');
if isfinite(bestCost)
    fprintf('Saved: qubo_profile_best_schedule.csv, qubo_profile_best_hourly.csv\n');
end
fprintf('NOTE: More profile choices need not improve finite-budget annealing.\n');
fprintf('Means exclude infeasible runs; full cost/runtime trade-offs are not matched to LP.\n');
