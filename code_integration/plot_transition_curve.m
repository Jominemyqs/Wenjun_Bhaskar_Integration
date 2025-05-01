function plot_transition_curve()
    baseDir  = 'LineSweep_JBB0p05';
    saveFile = fullfile(baseDir, 'pts.mat');
    figFile  = fullfile(baseDir, 'transition_curve_summary.png');

    if ~isfile(saveFile)
        error('Cannot find pts.mat in %s', baseDir);
    end

    load(saveFile, 'pts');

    figure(20); clf;
    plot(pts(:,1), pts(:,2), 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    xlabel('J_{BB}', 'FontSize', 12);
    ylabel('J_{OO}', 'FontSize', 12);
    title('Transition Curve (based on current simulations)', 'FontSize', 14);
    axis tight; grid on;

    saveas(gcf, figFile);
    fprintf('Saved curve summary to %s\n', figFile);
end
