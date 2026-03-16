function configurePlots(recon,test_num,cells,cfl,extra_prim_var_plot,riemann_solver,CD_Term_Order,plot_folder)
% configurePlots - Configure plot formatting, axis limits, and save EPS files.
%
% Uses a table-driven plot configuration to avoid hard-coded lookup logic.

if nargin < 8 || isempty(plot_folder)
    plot_folder = '../thesis/FigsforCommittee/';
end

file_pre = makeFilePrefix(recon, test_num, cells, cfl, riemann_solver, CD_Term_Order);
config = plotConfig(test_num, CD_Term_Order);
Folder = plot_folder;

% Standard figure styling (shared across all plots)
baseStyle = struct(
    'FontSize', 55,
    'LineWidth', 10,
    'MarkerSize', 4,
    'LegendFontSize', 60);

for figIdx = 1:length(config.vars)
    varCfg = config.vars{figIdx};

    figure(figIdx);
    set(gcf, 'Position', get(0, 'Screensize'));
    ax = gca;
    ax.YAxis.FontSize = baseStyle.FontSize;
    ax.XAxis.FontSize = baseStyle.FontSize;
    ax.LineWidth = baseStyle.LineWidth;

    % Apply y-axis limits if configured
    if ~any(isnan(varCfg.ylims(:)))
        ax.YAxis.Limits = varCfg.ylims(test_num,:);
    end

    % Legend
    legend('FontSize', baseStyle.LegendFontSize, 'Interpreter', 'Latex');
    legend('boxoff');
    leg = legend();
    leg.ItemTokenSize = [80, 25];
    leg.Position = config.locs(figIdx, :);

    % Labels
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 100, 'FontWeight', 'bold');
    ylabel(varCfg.ylabel, 'Interpreter', 'latex', 'FontSize', 100, 'FontWeight', 'bold');

    grid off;

    filename = fullfile(Folder, [file_pre, varCfg.suffix, '.eps']);
    saveas(gcf, filename, 'epsc');

    % Stop after optional vars if not requested
    if figIdx == 4 && ~extra_prim_var_plot
        break;
    end
end
end
