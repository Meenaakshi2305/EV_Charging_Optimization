%% Stage 4: Quantum-inspired EV charging using a genuine candidate-profile QUBO
% Requires MATLAB R2026a base MATLAB. No quantum hardware or extra packages.
% Place beside the 50-EV CSV, INR tariff CSV, and validated reference CSV.
% Model: choose exactly ONE complete charging profile per EV (binary choice).
% This is a restricted/discretized optimization, unlike the continuous LP.
clear; clc; close all;
rng(42,'twister'); % repeatable simulated annealing
EV = readtable('ev_charging_optimization_50_evs_synthetic.csv','TextType','string');
Tariff = readtable('ev_charging_hourly_tariffs_24h_inr.csv');
Ref = readtable('ev_charging_feasible_schedule_validation.csv','TextType','string');
N = height(EV); H = 24; cap = 250; dt = 1; tol = 1e-7;
assert(N==50 && height(Tariff)==24 && height(Ref)==N);
assert(isequal(Tariff.hour,(0:23)'), 'Tariffs must list hours 0:23.');
assert(isequal(string(EV.EV_ID),string(Ref.EV_ID)), 'Reference EV_ID row order mismatch.');
assert(all(EV.station_capacity_kw==cap), 'Unexpected station limit.');
price = Tariff.tariff_inr_per_kwh(:);
refPower = Ref{:,2:25};
assert(size(refPower,2)==H);
assert(all(abs(sum(refPower,2)*dt-EV.required_energy_kwh)<1e-5), 'Reference energy mismatch.');
assert(all(sum(refPower,1)<=cap+tol), 'Reference station violation.');

%% Build alternative charging profiles per EV (each already meets full EV demand).
% REFERENCE is always included. Other profiles ignore shared capacity initially.
% IMPORTANT: This is a FEASIBILITY-SAFEGUARDED DEMONSTRATION, not an independent
% competition with LP: the reference profile is known to cost INR 45190.05
% and is used to initialize the search. It therefore cannot establish
% whether annealing independently discovers the optimum.
% Candidate choices: reference, earliest, latest, cheapest, and 2 randomized.
% All charging amounts are in 0.1 kWh increments, matching input demands.
profiles = cell(N,1);
for i=1:N
    available = find((0:H-1)>=EV.arrival_hour(i) & (0:H-1)<EV.departure_hour(i));
    demand = round(10*EV.required_energy_kwh(i)); % integer units of 0.1 kWh
    maxUnits = floor(10*EV.max_charging_power_kw(i)+tol);
    assert(demand<=numel(available)*maxUnits,'Demand cannot fit the EV window.');
    P = round(refPower(i,:)*10)/10;
    options = P; % matrix rows: candidate profiles
    for k=1:5
        if k==1
            ordered = available; % earliest
        elseif k==2
            ordered = fliplr(available); % latest
        elseif k==3
            [~,idx] = sortrows([price(available), available(:)],[1 2]);
            ordered = available(idx); % cheapest hours
        else
            ordered = available(randperm(numel(available))); % diverse profiles
        end
        pUnits = zeros(1,H); remaining = demand;
        for h=ordered(:)'
            take = min(remaining,maxUnits);
            pUnits(h) = take;
            remaining = remaining-take;
            if remaining==0, break; end
        end
        assert(remaining==0);
        candidate = pUnits/10;
        if ~any(all(abs(options-candidate)<tol,2))
            options(end+1,:) = candidate; %#ok<SAGROW>
        end
    end
    assert(all(abs(sum(options,2)*dt-EV.required_energy_kwh(i))<tol));
    profiles{i}=options;
end

%% Explicit QUBO mathematical formulation
% Binary y(i,k): 1 if profile k is chosen for EV i.
% Binary s(h,b): slack (unused station capacity) at hour h.
% Slack in 0.1-kW increments via weights [1 2 ... 1024 453] * 0.1 kW.
% These weights represent every integer from 0 through 2500 (0-250 kW).
% E_QUBO = sum(profile_cost*y)
%        + A*sum_i(sum_k y(i,k)-1)^2
%        + B*sum_h(sum_i,k profile_power(i,k,h)*y(i,k)
%                   + sum_b slack_weight(b)*s(h,b) - 250)^2.
% All expressions are degree <= 2 in binary variables: a genuine QUBO.
% Simulated annealing uses group-wise binary swaps (one-hot preserved), and
% chooses slack bits analytically to minimize the QUBO for each proposed y.
% This is a restricted QUBO annealer, NOT a general-purpose QUBO solver.
slackWeights = [2.^(0:10), 453]/10;
assert(abs(sum(slackWeights)-cap)<tol);
A = 1e6; B = 1000; % INR-equivalent penalty weights; tune and report
choices = ones(N,1); % reference option for each EV; feasible starting state
loadKW = sum(refPower,1);
costs = cell(N,1);
for i=1:N, costs{i} = profiles{i}*price*dt; end
currentCost = sum(cellfun(@(c)c(1),costs));
currentPenalty = B*sum(max(0,loadKW-cap).^2);
bestChoices = choices; bestCost = currentCost; bestLoad = loadKW;

%% QUBO simulated annealing (group flips with optimized binary slack)
nRuns=8; stepsPerRun=25000;
bestFeasibleCost = bestCost; bestFeasibleChoices = choices;
bestFeasibleLoad = loadKW;
tic;
for run=1:nRuns
    % Keep a fresh feasible seed; improve best found across runs.
    if run==1
        choices = ones(N,1);
    else
        choices = bestFeasibleChoices;
    end
    loadKW = zeros(1,H); currentCost=0;
    for i=1:N
        loadKW=loadKW+profiles{i}(choices(i),:);
        currentCost=currentCost+costs{i}(choices(i));
    end
    currentPenalty = B*sum(max(0,loadKW-cap).^2);
    for t=1:stepsPerRun
        Tstart=4000; Tend=0.05;
        temperature = Tstart*(Tend/Tstart)^((t-1)/(stepsPerRun-1));
        i=randi(N);
        m=size(profiles{i},1);
        if m<2, continue; end
        proposed = randi(m-1);
        if proposed>=choices(i), proposed=proposed+1; end
        newLoad = loadKW - profiles{i}(choices(i),:) + profiles{i}(proposed,:);
        newPenalty = B*sum(max(0,newLoad-cap).^2);
        newCost = currentCost-costs{i}(choices(i))+costs{i}(proposed);
        dE = (newCost+newPenalty)-(currentCost+currentPenalty);
        if dE<=0 || rand < exp(-dE/temperature)
            choices(i)=proposed; loadKW=newLoad;
            currentCost=newCost; currentPenalty=newPenalty;
        end
        if all(loadKW<=cap+tol) && currentCost<bestFeasibleCost-1e-7
            bestFeasibleCost=currentCost;
            bestFeasibleChoices=choices;
            bestFeasibleLoad=loadKW;
        end
    end
    fprintf('Run %d/%d complete | best feasible INR %.2f\n',run,nRuns,bestFeasibleCost);
end
runtime=toc;

%% Verify final output against original EV requirements
chosenPower=zeros(N,H);
for i=1:N, chosenPower(i,:)=profiles{i}(bestFeasibleChoices(i),:); end
energy=sum(chosenPower,2)*dt;
station=sum(chosenPower,1);
evCosts=(chosenPower*price)*dt;
assert(all(abs(energy-EV.required_energy_kwh)<1e-6),'Energy validation failed.');
assert(all(station<=cap+1e-6),'Station capacity violated.');
for i=1:N
    outside=(0:H-1)<EV.arrival_hour(i) | (0:H-1)>=EV.departure_hour(i);
    assert(all(abs(chosenPower(i,outside))<tol));
    assert(all(chosenPower(i,:)<=EV.max_charging_power_kw(i)+tol));
end

% Verify that a binary slack encoding exists for all station hours.
% Each station load is a multiple of 0.1 kW; all unused capacities 0..250
% are representable by the listed slack bits.
slackBits=zeros(H,numel(slackWeights));
weights10=round(slackWeights*10);
for h=1:H
    remaining=round((cap-station(h))*10);
    for b=numel(weights10):-1:1
        if weights10(b)<=remaining
            slackBits(h,b)=1; remaining=remaining-weights10(b);
        end
    end
    assert(remaining==0,'Binary slack encoding failed.');
end
slackResidual = station(:) + slackBits*slackWeights(:) - cap;

assert(all(abs(slackResidual)<1e-6), ...
    'QUBO capacity-slack equality failed.');

% Verify actual QUBO energy: one-hot penalty zero and capacity penalty zero.
chosenCounts=zeros(N,1);
for i=1:N
    chosenCounts(i)=sum((1:size(profiles{i},1))==bestFeasibleChoices(i));
end
onehotPenalty=A*sum((chosenCounts-1).^2);
capacityPenalty=B*sum(slackResidual.^2);
quboEnergy=sum(evCosts)+onehotPenalty+capacityPenalty;

fprintf('\n--- QUBO QUANTUM-INSPIRED RESULTS ---\n');
fprintf('Algorithm: profile-choice QUBO with annealing (restricted model)\n');
fprintf('EVs: %d, candidate profiles: %d, slack bits: %d\n', N, sum(cellfun(@(p)size(p,1),profiles)), H*numel(slackWeights));
fprintf('Total energy delivered: %.1f kWh\n',sum(energy));
fprintf('Total electricity cost: INR %.2f\n',sum(evCosts));
fprintf('Peak station demand: %.2f kW\n',max(station));
fprintf('Capacity violations: %d\n',sum(station>cap+1e-6));
fprintf('Max unmet EV energy: %.6f kWh\n',max(abs(energy-EV.required_energy_kwh)));
fprintf('QUBO energy (including penalties): %.2f\n',quboEnergy);
fprintf('Annealing runtime: %.3f seconds\n',runtime);
fprintf('Continuous LP optimum: INR 45190.05 (from validated LP run)\n');
fprintf('Cost difference vs LP benchmark: %.2f INR\n',sum(evCosts)-45190.05);
fprintf('CAUTION: Feasible seed includes the supplied reference schedule.\n');
fprintf('This demonstration DOES NOT independently benchmark QUBO vs LP.\n');
fprintf('Status: FEASIBLE (not proven globally optimal for QUBO).\n');

EVResults=table(EV.EV_ID,EV.required_energy_kwh,energy,evCosts,bestFeasibleChoices, ...
 'VariableNames',{'EV_ID','Required_kWh','Delivered_kWh','Cost_INR','SelectedProfile'});
writetable(EVResults,'qubo_ev_results.csv');
HourlyResults=table((0:23)',station',repmat(cap,H,1),price, ...
 'VariableNames',{'hour','Station_Demand_kW','Station_Capacity_kW','Tariff_INR_per_kWh'});
writetable(HourlyResults,'qubo_hourly_results.csv');
Schedule=array2table(chosenPower,'VariableNames',compose('hour_%02d_kw',0:23));
Schedule=addvars(Schedule,EV.EV_ID,'Before',1,'NewVariableNames','EV_ID');
writetable(Schedule,'qubo_schedule.csv');
figure('Name','QUBO charging schedule');
plot(0:23,station,'-o','LineWidth',1.5); hold on;
yline(cap,'r--','250 kW station limit','LineWidth',1.5);
xlabel('Hour of day'); ylabel('Aggregate station demand (kW)');
title('Quantum-inspired QUBO: feasible charging demand'); grid on; xlim([0 23]);
