figure('Color','w');

titles = {'Uncontrolled Charging', ...
    'Classical LP', ...
    'Independent QUBO'};

schedules = {P_uncontrolled, P_lp, P_qubo};

for k = 1:3
    subplot(1,3,k);

    imagesc(0:23, 1:50, schedules{k});

    title(titles{k});
    xlabel('Hour of day');
    ylabel('EV Number');

    colormap(parula);
    colorbar;
    caxis([0, maxPower]);

    set(gca,'FontSize',9);
end

sgtitle('EV Charging Schedule Comparison');

exportgraphics(gcf, ...
    'charging_schedule_comparison.png', ...
    'Resolution',300);