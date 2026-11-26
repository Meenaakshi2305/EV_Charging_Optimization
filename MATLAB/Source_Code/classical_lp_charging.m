%% Stage 3 - Classical EV charging optimization using Linear Programming
% Place this script beside the EV and INR tariff CSV files.
clear; clc; close all;

EV = readtable('ev_charging_optimization_50_evs_synthetic.csv', 'TextType','string');
Tariff = readtable('ev_charging_hourly_tariffs_24h_inr.csv');

nEV = height(EV);
nH = 24;
slotHours = 1;
capacity = 250;
assert(nEV == 50 && height(Tariff) == nH, 'Expected 50 EVs and 24 tariff slots.');
assert(isequal(Tariff.hour, (0:23)'), 'Tariff hours must be 0 through 23 in order.');
assert(all(EV.station_capacity_kw == capacity), 'Station capacity is not consistently 250 kW.');
assert(all(EV.arrival_hour >= 0 & EV.departure_hour <= nH & EV.arrival_hour < EV.departure_hour), ...
    'Invalid EV arrival/departure data.');
assert(all(EV.required_energy_kwh > 0 & EV.max_charging_power_kw > 0), 'Invalid EV energy/power.');

% Variables x(i,h) = average power in kW supplied to EV i during hour h.
% MATLAB uses column-major order: x(:) stacks all EVs in hour 0, then hour 1, etc.
upperPower = zeros(nEV,nH);
for i = 1:nEV
    allowedHours = (0:nH-1) >= EV.arrival_hour(i) & ...
                   (0:nH-1) < EV.departure_hour(i);
    upperPower(i,allowedHours) = EV.max_charging_power_kw(i);
end

% Minimize sum of [power (kW) * duration (h) * tariff (INR/kWh)].
costCoefficient = repmat(Tariff.tariff_inr_per_kwh(:)' * slotHours, nEV, 1);
f = costCoefficient(:);

% Each EV receives exactly its required energy.
Aeq = kron(ones(1,nH), speye(nEV)) * slotHours;
beq = EV.required_energy_kwh;

% Every hour's combined power stays <= station capacity.
A = kron(speye(nH), ones(1,nEV));
b = capacity * ones(nH,1);

% All x are nonnegative; upper bounds enforce EV-specific max power and windows.
lb = zeros(nEV*nH,1);
ub = upperPower(:);

options = optimoptions('linprog', 'Display','none');
[xVector, totalCost, exitflag, output] = linprog(f,A,b,Aeq,beq,lb,ub,options);
assert(exitflag == 1, 'LP did not report an optimal solution: %s', output.message);

powerSchedule = reshape(xVector,nEV,nH);
energyDelivered = sum(powerSchedule,2) * slotHours;
hourlyDemand = sum(powerSchedule,1)';
individualCost = sum(powerSchedule .* costCoefficient,2);
tol = 1e-5;
assert(all(abs(energyDelivered - EV.required_energy_kwh) < tol), 'Some EV demands were not met.');
assert(all(hourlyDemand <= capacity + tol), 'Station capacity exceeded.');
assert(all(powerSchedule(:) >= -tol & powerSchedule(:) <= ub + tol), 'Power/window violation.');
assert(abs(sum(individualCost)-totalCost) < 1e-4, 'Cost mismatch.');

% Detailed outputs for report and comparison with other methods.
EVResults = table(EV.EV_ID, EV.required_energy_kwh, energyDelivered, individualCost, ...
    'VariableNames', {'EV_ID','Required_kWh','Delivered_kWh','Cost_INR'});
HourlyResults = table((0:nH-1)', Tariff.tariff_inr_per_kwh, hourlyDemand, ...
    capacity*ones(nH,1), 'VariableNames', ...
    {'Hour','Tariff_INR_per_kWh','StationDemand_kW','StationCapacity_kW'});
ScheduleResults = array2table(powerSchedule, ...
    'VariableNames', cellstr(compose('Hour_%02d',0:nH-1)));
ScheduleResults = addvars(ScheduleResults, EV.EV_ID, 'Before',1,'NewVariableNames','EV_ID');
writetable(EVResults,'classical_lp_ev_results.csv');
writetable(HourlyResults,'classical_lp_hourly_results.csv');
writetable(ScheduleResults,'classical_lp_schedule.csv');

fprintf('\n--- CLASSICAL LP OPTIMIZATION ---\n');
fprintf('Solver status: OPTIMAL\n');
fprintf('Total energy delivered: %.1f kWh\n',sum(energyDelivered));
fprintf('Total charging cost: INR %.2f\n',totalCost);
fprintf('Peak station demand: %.2f kW\n',max(hourlyDemand));
fprintf('Hours exceeding station capacity: %d\n',sum(hourlyDemand>capacity+tol));
fprintf('Maximum unmet EV energy: %.6f kWh\n',max(abs(energyDelivered-EV.required_energy_kwh)));
fprintf('Solver iterations: %d\n',output.iterations);

figure;
bar(0:nH-1,hourlyDemand, 'FaceColor',[0.2 0.5 0.8]); hold on;
yline(capacity,'r--','250 kW station limit','LineWidth',1.5);
xlabel('Hour of day (0 to 23)'); ylabel('Station charging demand (kW)');
title('Classical LP: optimized hourly charging demand');
xticks(0:23); grid on;
