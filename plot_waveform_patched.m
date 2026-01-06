function plot_waveform_patched(ax, X, Y, WF, cardiacCycle, clr, PLW)

patch(ax, X, Y, clr, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold(ax, 'on');
plot(ax, cardiacCycle, WF, 'LineWidth', PLW, 'Color', clr);
hold(ax, 'off');

end