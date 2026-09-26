%% Replot the hardware-consistent columns of the 10,000-trial pilot study
% A fixed reference phase at a representative receiving unit is the scalar
% restriction of the fixed spatial code. Time-varying random/four-phase columns
% from the source study are intentionally excluded.

clc; clear; close all;

script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');

T = readtable(fullfile(result_dir, 'pilot_length_source_summary.csv'));

names = {'No reference', 'Constant'};
colors = [.35 .35 .35; 0 .447 .741];

% Set defaults before creating figures.
set(groot, ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultTextFontName', 'Times New Roman', ...
    'defaultAxesFontSize', 11, ...
    'defaultLineLineWidth', 1.5, ...
    'defaultLineMarkerSize', 6);

% Output names for the two separate figures.
file_names = { ...
    'fig_spatial_pilot_recovery_probability', ...
    'fig_spatial_pilot_channel_nmse'};

for panel = 1:2

    % Create an independent figure for each panel.
    fig = figure('Color', 'w', 'Position', [100 100 450 340]);
    ax = axes('Parent', fig);
    hold(ax, 'on');

    h = gobjects(2, 1);

    for k = 1:2
        idx = strcmp(T.ReferenceDesign, names{k});
        Tk = T(idx, :);

        % Sort pilot lengths to ensure correct line connections.
        Tk = sortrows(Tk, 'PilotLength');

        if panel == 1
            h(k) = plot(ax, ...
                Tk.PilotLength, Tk.NoiselessRecoveryProbability, ...
                'o-', 'Color', colors(k, :));
        else
            h(k) = plot(ax, ...
                Tk.PilotLength, Tk.ChannelNMSEdB, ...
                'o-', 'Color', colors(k, :));
        end
    end

    if panel == 1
        ylabel(ax, 'Recovery probability');
        ylim(ax, [0 1.03]);
        % title(ax, 'Noiseless multistart recovery');
    else
        ylabel(ax, 'Channel NMSE (dB)');
        % title(ax, 'Noisy recovery at 12 dB');
    end

    xline(ax, 2, ':', 'HandleVisibility', 'off');
    xline(ax, 4, ':', 'HandleVisibility', 'off');
    xline(ax, 6, ':', 'HandleVisibility', 'off');

    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, 'Pilot symbols, T_p');
    xlim(ax, [1 48]);

    % Add a legend to each figure, using only the two plotted curves.
    legend(ax, h, names, 'Location', 'best', 'FontSize', 9);

    % Save each figure separately.
    exportgraphics(fig, ...
        fullfile(result_dir, [file_names{panel}, '.eps']), ...
        'ContentType', 'vector');

    exportgraphics(fig, ...
        fullfile(result_dir, [file_names{panel}, '.png']), ...
        'Resolution', 300);

    savefig(fig, fullfile(result_dir, [file_names{panel}, '.fig']));

    close(fig); % 注释掉此行，可在运行结束后保留两个图窗。
end

% Export the selected data.
idx = strcmp(T.ReferenceDesign, 'No reference') | ...
      strcmp(T.ReferenceDesign, 'Constant');

writetable(T(idx, :), ...
    fullfile(result_dir, 'spatial_pilot_identifiability_summary.csv'));