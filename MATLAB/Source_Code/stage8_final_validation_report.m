%% STAGE 8 - AUDITED COMPARISON (MATLAB R2026a)
% Reads existing algorithm-produced CSVs. Never inserts the QUBO optimum.
% If a required result is missing, re-run that stage first.
clear; clc; close all;
tol = 1e-5; dt = 1;
EV = readtable('ev_charging_optimization_50_evs_synthetic.csv','TextType','string');
T = readtable('ev_charging_hourly_tariffs_24h_inr.csv');
assert(height(EV)==50 && height(T)==24,'Unexpected dataset dimensions.');
assert(isequal(T.hour(:),(0:23)'),'Tariff hours must be 0 through 23.');
assert(numel(unique(EV.EV_ID))==50,'EV_ID must be unique.');
assert(all(EV.station_capacity_kw==250),'Station capacity is not uniform at 250 kW.');
price = T.tariff_inr_per_kwh(:); cap = 250;
assert(all(isfinite(price) & price>=0),'Invalid hourly tariffs.');
assert(all(EV.required_energy_kwh>0 & EV.max_charging_power_kw>0),'Invalid EV energy/power.');
assert(all(EV.arrival_hour>=0 & EV.departure_hour<=24 & EV.arrival_hour<EV.departure_hour),'Invalid EV windows.');

methodNames = {'Uncontrolled';'Classical LP';'Independent QUBO'};
filenames = {'uncontrolled_schedule.csv'; 'classical_lp_schedule.csv'; 'qubo_profile_best_schedule.csv'};
N=50; H=24;
costs=zeros(3,1); peaks=zeros(3,1); excessHours=zeros(3,1);
maxExcess=zeros(3,1); maxUnmet=zeros(3,1); stationFeasible=false(3,1);
demands=zeros(H,3); charges=cell(3,1);
for m=1:3
    assert(isfile(filenames{m}), ['Missing ',filenames{m},'. Run its original MATLAB script first.']);
    S=readtable(filenames{m},'TextType','string');
    assert(height(S)==N && width(S)==H+1, ['Unexpected schedule shape: ',filenames{m}]);
    assert(isequal(string(S.EV_ID),string(EV.EV_ID)), ['EV order/IDs mismatch: ',filenames{m}]);
    P=S{:,2:end};
    assert(isnumeric(P) && all(isfinite(P(:))),'Non-numeric or missing power values.');
    assert(all(P(:)>=-tol),'Negative charging power.');
    assert(all(all(P <= EV.max_charging_power_kw+tol)),'Per-EV maximum power exceeded.');
    for i=1:N
        forbidden=(0:H-1)<EV.arrival_hour(i) | (0:H-1)>=EV.departure_hour(i);
        assert(all(abs(P(i,forbidden))<tol),['Charging outside EV window in ',filenames{m}]);
    end
    delivered=sum(P,2)*dt;
    maxUnmet(m)=max(abs(delivered-EV.required_energy_kwh));
    assert(maxUnmet(m)<tol,['Energy requirement not met in ',filenames{m}]);
    demand=sum(P,1)' ;
    demands(:,m)=demand; charges{m}=P;
    costs(m)=sum(demand.*price)*dt;
    peaks(m)=max(demand);
    excessHours(m)=sum(demand>cap+tol);
    maxExcess(m)=max(0,max(demand-cap));
    stationFeasible(m)=excessHours(m)==0;
end

% These are regressions against previously observed MATLAB logs, NOT
% substituted results. A mismatch stops the report instead of hiding it.
assert(abs(costs(1)-46450.80)<0.02,'Uncontrolled cost differs from verified run.');
assert(abs(costs(2)-45190.05)<0.02,'Classical LP cost differs from verified run.');
assert(abs(costs(3)-45996.20)<0.02,'Stage 7 QUBO cost differs. Verify selected output CSV.');
assert(excessHours(1)==8 && all(excessHours(2:3)==0),'Unexpected capacity violations.');
assert(abs(peaks(1)-571.3)<0.02,'Unexpected uncontrolled peak.');
assert(all(abs(sum(charges{1},2)-EV.required_energy_kwh)<tol));

% Validate sensitivity experiment summaries directly from the run data.
assert(isfile('qubo_budget_all_runs.csv') && isfile('qubo_profile_all_runs.csv'), ...
    'Missing experiment CSV files. Run Stages 6 and 7 first.');
B=readtable('qubo_budget_all_runs.csv');
R=readtable('qubo_profile_all_runs.csv');
assert(height(B)==80 && height(R)==60,'Unexpected experiment run counts.');
assert(all(ismember([5000 15000 30000 60000],unique(B.BudgetSteps))) && ...
       all(ismember([12 24 48],unique(R.GroupLimit))),'Unexpected experiment groups.');

if isfile('qubo_budget_summary.csv')
    budgetSummary=readtable('qubo_budget_summary.csv');
else
    error('Missing qubo_budget_summary.csv. Re-run Stage 6.');
end
if isfile('qubo_profile_summary.csv')
    profileSummary=readtable('qubo_profile_summary.csv');
else
    error('Missing qubo_profile_summary.csv. Re-run Stage 7.');
end
for k=1:height(budgetSummary)
    g=B(B.BudgetSteps==budgetSummary.Budget(k),:);
    good=g.FeasibleCost(g.Feasible==1);
    assert(height(g)==20 && numel(good)==budgetSummary.FeasibleRuns(k),'Stage 6 feasibility mismatch.');
    assert(abs(min(good)-budgetSummary.BestFeasibleCost_INR(k))<.02,'Stage 6 best cost mismatch.');
end
for k=1:height(profileSummary)
    g=R(R.GroupLimit==profileSummary.CandidateLimit(k),:);
    good=g.FeasibleCost_INR(g.Feasible==1);
    assert(height(g)==20 && numel(good)==profileSummary.FeasibleRuns(k),'Stage 7 feasibility mismatch.');
    assert(abs(min(good)-profileSummary.BestFeasibleCost_INR(k))<.02,'Stage 7 best cost mismatch.');
end
assert(abs(min(R.FeasibleCost_INR(R.Feasible==1))-costs(3))<.02, ...
    'Saved Stage 7 schedule does not match best Stage 7 run.');

Comparison=table(string(methodNames),costs,repmat(sum(EV.required_energy_kwh),3,1),peaks,excessHours, ...
    maxExcess,maxUnmet,stationFeasible,'VariableNames', ...
    {'Method','ElectricityCost_INR','DeliveredEnergy_kWh','PeakDemand_kW', ...
     'CapacityViolationHours','MaximumOverload_kW','MaximumEnergyError_kWh','StationFeasible'});
writetable(Comparison,'stage8_method_comparison.csv');
Hourly=table((0:23)',price,demands(:,1),demands(:,2),demands(:,3),repmat(cap,H,1), ...
    'VariableNames',{'Hour','Tariff_INR_per_kWh','Uncontrolled_kW','ClassicalLP_kW','IndependentQUBO_kW','Capacity_kW'});
writetable(Hourly,'stage8_hourly_comparison.csv');
writetable(budgetSummary,'stage8_annealing_budget_summary.csv');
writetable(profileSummary,'stage8_profile_sensitivity_summary.csv');

figure('Name','Stage 8: cost comparison','Color','w');
bar(costs); set(gca,'XTick',1:3,'XTickLabel',methodNames);
ylabel('Charging cost (INR)'); title('Cost comparison: uncontrolled schedule is infeasible'); grid on;
exportgraphics(gcf,'stage8_cost_comparison.png','Resolution',200);

figure('Name','Stage 8: hourly power','Color','w');
plot(0:23,demands(:,1),'-o','LineWidth',1.4); hold on;
plot(0:23,demands(:,2),'-s','LineWidth',1.4);
plot(0:23,demands(:,3),'-d','LineWidth',1.4);
yline(cap,'--','Station capacity (250 kW)','LineWidth',1.3);
xlabel('Hour (start of slot)'); ylabel('Station charging demand (kW)');
title('Hourly station demand'); legend('Uncontrolled','Classical LP','Independent QUBO','Capacity','Location','best');
xticks(0:23); grid on;
exportgraphics(gcf,'stage8_hourly_demand.png','Resolution',200);

figure('Name','Stage 8: sensitivity','Color','w');
plot(budgetSummary.Budget,budgetSummary.MeanFeasibleCost_INR,'-o','LineWidth',1.5); hold on;
yline(costs(2),'--','Continuous LP optimum');
xlabel('Annealing steps per run'); ylabel('Mean feasible QUBO cost (INR)');
title('Annealing budget sensitivity (feasible runs only)'); grid on;
exportgraphics(gcf,'stage8_annealing_sensitivity.png','Resolution',200);

figure('Name','Stage 8: profile diversity','Color','w');
plot(profileSummary.UniqueProfiles,profileSummary.MeanFeasibleCost_INR,'-o','LineWidth',1.5); hold on;
plot(profileSummary.UniqueProfiles,profileSummary.BestFeasibleCost_INR,'-s','LineWidth',1.5);
yline(costs(2),'--','Continuous LP optimum');
xlabel('Number of unique candidate profiles'); ylabel('QUBO charging cost (INR)');
title('Candidate profile sensitivity'); legend('Mean feasible','Best feasible','LP optimum','Location','best'); grid on;
exportgraphics(gcf,'stage8_profile_sensitivity.png','Resolution',200);

fprintf('\n--- STAGE 8 INDEPENDENT CSV AUDIT ---\n'); disp(Comparison);
fprintf('Classical LP savings relative to infeasible uncontrolled baseline: INR %.2f (%.3f%%)\n', ...
    costs(1)-costs(2),100*(1-costs(2)/costs(1)));
fprintf('Best independent QUBO gap to LP: INR %.2f (%.3f%%)\n', ...
    costs(3)-costs(2),100*(costs(3)/costs(2)-1));
fprintf('Stage 6: %d/%d feasible; Stage 7: %d/%d feasible.\n', ...
    sum(B.Feasible==1),height(B),sum(R.Feasible==1),height(R));
fprintf('WARNING: LP and QUBO are different decision spaces (continuous power vs selected profiles).\n');
fprintf('WARNING: This script does not compare end-to-end LP and QUBO runtimes.\n');
fprintf('WARNING: Peak demand is capacity-constrained, not explicitly minimized.\n');
fprintf('PASS: All input/output schedule audits completed.\n');
