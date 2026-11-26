%% STAGE 2: Uncontrolled EV charging baseline
% One-hour slots; arrival inclusive, departure exclusive.
% EVs start at arrival and charge as fast as possible.
% Station capacity is NOT enforced here: violations are reported.
clear; clc; close all;

%% 1. Read the exact same input files used in validation
EV = readtable('ev_charging_optimization_50_evs_synthetic.csv', 'TextType','string');
Tariff = readtable('ev_charging_hourly_tariffs_24h_inr.csv');

%% 2. Confirm essential inputs
assert(height(EV)==50, 'Expected exactly 50 EVs.');
assert(height(Tariff)==24 && isequal(Tariff.hour(:),(0:23)'), ...
    'Tariff hours must be 0 through 23 in order.');
assert(all(abs(EV.station_capacity_kw-250)<1e-9), ...
    'Expected shared station capacity of 250 kW.');
assert(numel(unique(EV.EV_ID))==height(EV), 'EV IDs are not unique.');
assert(all(EV.required_energy_kwh>0 & EV.max_charging_power_kw>0), ...
    'Required energy and maximum power must be positive.');
assert(all(EV.arrival_hour>=0 & EV.departure_hour<=24 & ...
    EV.arrival_hour<EV.departure_hour & ...
    EV.arrival_hour==floor(EV.arrival_hour) & ...
    EV.departure_hour==floor(EV.departure_hour)), 'Invalid EV time windows.');
assert(all(EV.required_energy_kwh <= EV.max_charging_power_kw .* ...
    (EV.departure_hour-EV.arrival_hour)+1e-8), ...
    'Some EVs cannot meet their demand in their availability window.');
assert(all(isfinite(Tariff.tariff_inr_per_kwh) & Tariff.tariff_inr_per_kwh>=0), ...
    'Invalid tariffs.');

%% 3. Make the charging schedule (rows = EVs; columns = hours 0..23)
nEV = height(EV);
hours = 0:23;
slotHours = 1;
powerSchedule = zeros(nEV,24);  % kW (each slot lasts one hour)
for i = 1:nEV
    energyLeft = EV.required_energy_kwh(i); % kWh
    for h = EV.arrival_hour(i):(EV.departure_hour(i)-1)
        % Charge immediately at the highest allowed rate, until full.
        p = min(EV.max_charging_power_kw(i), energyLeft/slotHours);
        powerSchedule(i,h+1) = p; % MATLAB column 1 corresponds to hour 0
        energyLeft = energyLeft - p*slotHours;
        if energyLeft <= 1e-9
            break;
        end
    end
end

%% 4. Check delivered energy and compute each EV's cost
prices = Tariff.tariff_inr_per_kwh(:); % 24 x 1, INR/kWh
energyDelivered = sum(powerSchedule,2)*slotHours; % 50 x 1, kWh
costPerEV = (powerSchedule*prices)*slotHours;  % 50 x 1, INR
unmetEnergy = max(0,EV.required_energy_kwh-energyDelivered);
assert(all(unmetEnergy<1e-6), 'Some EVs did not receive required energy.');

%% 5. Station demand and violations (do NOT silently enforce limit)
stationLimit = 250; % kW
hourlyDemand = sum(powerSchedule,1); % 1 x 24, kW
violating = hourlyDemand > stationLimit+1e-8;
hoursAboveLimit = sum(violating);
peakDemand = max(hourlyDemand);
maxExcess = max(0,peakDemand-stationLimit);
totalCost = sum(costPerEV);

%% 6. Results tables and console summary
results = table(EV.EV_ID,EV.required_energy_kwh,energyDelivered, ...
    costPerEV,unmetEnergy, ...
    'VariableNames',{'EV_ID','Required_kWh','Delivered_kWh','Cost_INR','Unmet_kWh'});
hourlyResults = table(hours(:),hourlyDemand(:), ...
    repmat(stationLimit,24,1),violating(:),prices, ...
    'VariableNames',{'Hour','Demand_kW','Limit_kW','ExceedsLimit','Tariff_INR_per_kWh'});

disp('--- INDIVIDUAL EV CHARGING RESULTS ---');
disp(results);
fprintf('\n--- UNCONTROLLED BASELINE SUMMARY ---\n');
fprintf('EVs: %d\n',nEV);
fprintf('Total energy delivered: %.1f kWh\n',sum(energyDelivered));
fprintf('Total cost (all EVs): INR %.2f\n',totalCost);
fprintf('Peak station demand: %.2f kW\n',peakDemand);
fprintf('Station capacity: %.2f kW\n',stationLimit);
fprintf('Hours exceeding capacity: %d\n',hoursAboveLimit);
fprintf('Maximum capacity excess: %.2f kW\n',maxExcess);
if hoursAboveLimit > 0
    fprintf('STATUS: INFEASIBLE with the shared 250 kW station limit.\n');
else
    fprintf('STATUS: Within shared station capacity.\n');
end

%% 7. Save outputs in separate files
writetable(results,'uncontrolled_ev_results.csv');
writetable(hourlyResults,'uncontrolled_hourly_results.csv');
scheduleTable = array2table(powerSchedule, ...
    'VariableNames',compose('hour_%02d',hours));
scheduleTable = addvars(scheduleTable,EV.EV_ID,'Before',1, ...
    'NewVariableNames','EV_ID');
writetable(scheduleTable,'uncontrolled_schedule.csv');

%% 8. Plots
figure('Name','Uncontrolled charging: station demand');
bar(hours,hourlyDemand);
hold on;
yline(stationLimit,'r--','250 kW station limit','LineWidth',1.5);
hold off;
xlabel('Hour of day (0 to 23)');
ylabel('Aggregate charging power (kW)');
title('Uncontrolled EV charging: hourly station demand');
xticks(hours); grid on;

figure('Name','Uncontrolled charging: cost by EV');
bar(costPerEV);
xlabel('EV number (1 to 50)');
ylabel('Charging cost (INR)');
title('Uncontrolled EV charging: individual costs');
grid on;
