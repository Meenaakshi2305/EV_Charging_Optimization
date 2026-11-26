%% Quantum-Inspired EV Charging Optimization - Input Validation
% MATLAB R2026a. Place this script and the three CSVs in the same folder.
clear; clc;
projectFolder = fileparts(mfilename('fullpath'));
evFile = fullfile(projectFolder,'ev_charging_optimization_50_evs_synthetic.csv');
tariffFile = fullfile(projectFolder,'ev_charging_hourly_tariffs_24h_inr.csv');
scheduleFile = fullfile(projectFolder,'ev_charging_feasible_schedule_validation.csv');
assert(isfile(evFile) && isfile(tariffFile) && isfile(scheduleFile), ...
    'A required CSV is missing. Save all three CSVs beside this script.');

% Import data. The original uploaded EV dataset is not modified.
EV = readtable(evFile, 'TextType','string');
Tariff = readtable(tariffFile);
Schedule = readtable(scheduleFile, 'TextType','string');

% Verify the exact expected columns and row counts.
evColumns = ["EV_ID", "required_energy_kwh", "max_charging_power_kw", ...
    "arrival_hour", "departure_hour", "station_capacity_kw"];
assert(all(ismember(evColumns,string(EV.Properties.VariableNames))), ...
    'EV file has missing expected columns.');
assert(all(ismember(["hour","tariff_inr_per_kwh"], ...
    string(Tariff.Properties.VariableNames))), 'Tariff columns are incorrect.');
assert(height(EV)==50,'EV file must contain exactly 50 EVs.');
assert(height(Tariff)==24,'Tariff file must have 24 hourly prices.');
assert(numel(unique(EV.EV_ID))==50 && all(strlength(EV.EV_ID)>0), ...
    'EV_ID values must be distinct and nonempty.');
assert(~any(ismissing(EV),'all') && ~any(ismissing(Tariff),'all'), ...
    'Missing values found in EV or tariff data.');

% Check that the numeric inputs are finite and valid.
numbers = [EV.required_energy_kwh, EV.max_charging_power_kw, ...
    EV.arrival_hour, EV.departure_hour, EV.station_capacity_kw];
assert(all(isfinite(numbers),'all'),'EV numerical fields must be finite.');
assert(all(EV.required_energy_kwh>0) && all(EV.max_charging_power_kw>0), ...
    'Required energy and charging-power limits must be positive.');
assert(all(EV.arrival_hour>=0 & EV.departure_hour<=24 & ...
    EV.arrival_hour<EV.departure_hour), 'Invalid arrival/departure windows.');
assert(all(mod(EV.arrival_hour,1)==0 & mod(EV.departure_hour,1)==0), ...
    'Arrival and departure must be integer hours.');
assert(all(EV.station_capacity_kw==250),'Station capacity must be 250 kW.');
assert(all(EV.required_energy_kwh <= ...
    EV.max_charging_power_kw.*(EV.departure_hour-EV.arrival_hour)+1e-6), ...
    'At least one EV has an individually impossible energy request.');
assert(all(Tariff.hour==(0:23)'), 'Tariff hours must be 0, 1, ..., 23.');
assert(all(isfinite(Tariff.tariff_inr_per_kwh) & ...
    Tariff.tariff_inr_per_kwh>=0), 'Tariff entries must be finite and nonnegative.');

% Verify the provided reference schedule, which is NOT an optimal solution.
slotNames = "hour_" + compose('%02d',0:23) + "_kw";
assert(height(Schedule)==50 && all(ismember(["EV_ID",slotNames], ...
    string(Schedule.Properties.VariableNames))), 'Reference schedule layout is incorrect.');
assert(numel(unique(Schedule.EV_ID))==50 && ...
    all(ismember(EV.EV_ID,Schedule.EV_ID)), 'Schedule EV IDs do not match EV data.');
[matched, idx] = ismember(EV.EV_ID,Schedule.EV_ID);
assert(all(matched),'Unmatched schedule IDs.');
P = Schedule{idx,cellstr(slotNames)}; % Rows: EV order; columns: hours 0..23
assert(isnumeric(P) && all(isfinite(P),'all'),'Invalid power values in schedule.');
tolerance = 1e-5;
assert(all(P>=-tolerance,'all'),'Negative charging power detected.');
assert(all(P<=EV.max_charging_power_kw+tolerance,'all'), ...
    'Reference schedule exceeds an individual EV power limit.');
H = 0:23;
available = H>=EV.arrival_hour & H<EV.departure_hour;
assert(all(abs(P(~available))<=tolerance), ...
    'Reference schedule has charging outside an availability window.');
assert(all(abs(sum(P,2)-EV.required_energy_kwh)<=tolerance), ...
    'At least one EV has an energy deficit or excess.');
stationLoad = sum(P,1);
assert(all(stationLoad<=250+tolerance), ...
    'Reference schedule exceeds the station capacity.');
referenceCostINR = sum(stationLoad .* Tariff.tariff_inr_per_kwh');

fprintf('PASS: %d distinct EVs and %d tariff hours imported.\n',height(EV),height(Tariff));
fprintf('PASS: Total requested energy = %.1f kWh.\n',sum(EV.required_energy_kwh));
fprintf('PASS: Reference schedule meets energy, window, and power constraints.\n');
fprintf('Reference peak power: %.1f kW (capacity: 250 kW).\n',max(stationLoad));
fprintf('Reference schedule cost: INR %.2f (NOT an optimized cost).\n',referenceCostINR);

figure;
plot(Tariff.hour,Tariff.tariff_inr_per_kwh,'-o','LineWidth',1.5);
xlabel('Hour of day (0 to 23)'); ylabel('Tariff (INR/kWh)');
title('Illustrative 24-hour time-of-use electricity tariff');
grid on; xlim([0 23]); xticks(0:23);

