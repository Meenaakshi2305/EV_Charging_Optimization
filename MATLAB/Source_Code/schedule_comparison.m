%% CHARGING SCHEDULE HEATMAP COMPARISON - 50 EVs, 24 hours
% Save in the folder containing the three schedule CSV files.
clear; clc; close all;

files = {'uncontrolled_schedule.csv', ...
         'classical_lp_schedule.csv', ...
         'qubo_profile_best_schedule.csv'};
labels = {'Uncontrolled Charging', 'Classical LP', 'Independent QUBO'};
P = cell(1,3);
IDs = strings(50,1);

for k = 1:3
    assert(isfile(files{k}), 'Cannot find %s. Run the corresponding model first.', files{k});
    data = readtable(files{k}, 'TextType','string');
    assert(height(data)==50 && width(data)==25, ...
        '%s must contain 50 EV rows, one EV_ID column and 24 power columns.',files{k});
    assert(ismember('EV_ID',data.Properties.VariableNames), 'EV_ID column missing in %s.',files{k});
    if k == 1
        IDs = string(data.EV_ID);
    else
        assert(isequal(IDs,string(data.EV_ID)), 'EV order differs in %s.', files{k});
    end
    P{k} = data{:,2:25};
    assert(isnumeric(P{k}) && all(isfinite(P{k}(:))), ...
        'Charging powers must be numeric in %s.',files{k});
    fprintf('%s: %d EVs, %d hours, total %.1f kWh, peak %.1f kW\n', ...
        labels{k},size(P{k},1),size(P{k},2),sum(P{k}(:)),max(sum(P{k},1)));
end

maxPower = max(cellfun(@(x) max(x(:)),P));
fig = figure('Color','w','Position',[100 100 1600 680]);
t = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
for k = 1:3
    ax = nexttile(t);
    imagesc(ax,0:23,1:50,P{k});
    set(ax,'YDir','normal');
    caxis(ax,[0 maxPower]);
    colormap(ax,parula);
    title(ax,labels{k},'FontWeight','bold');
    xlabel(ax,'Hour of day');
    if k == 1
        ylabel(ax,'EV number (1-50)');
    end
    xticks(ax,0:3:23);
    yticks(ax,1:5:50);
    ax.FontSize = 10;
end
cb = colorbar('Location','eastoutside');
cb.Layout.Tile = 'east';
cb.Label.String = 'Charging power (kW)';
title(t,'EV Charging Schedule Comparison - 50 EVs, 24 Hours');
exportgraphics(fig,'charging_schedule_comparison.png','Resolution',300);
fprintf('Saved: charging_schedule_comparison.png\n');
