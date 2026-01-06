function plot_waveform_lines(ax, WF, cardiacCycle, clr, PLW)

plot(ax, cardiacCycle, WF.avg, 'LineWidth', PLW*2, 'Color', clr);
hold(ax, 'on');
for k = 1:size(WF.all, 2)
    plot(ax, cardiacCycle, WF.all(:, k), ...
        'LineWidth', PLW/2, 'Color', c1, ...
        'LineStyle', '--', ...
        'Marker', 'none');
end
ax.YColor = clr;

end