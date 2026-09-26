%% Replot stability-guided training as two independent figures
% This script loads the saved Monte Carlo results and changes only the
% presentation. No channel generation or Monte Carlo trial is rerun.
clc;
clear;
close all;

script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');
data = load(fullfile(result_dir, ...
    'stability_guided_training_baseline_results.mat'));

snr_dB = data.snr_dB;
method_names = data.method_names;
channel_nmse_dB = data.channel_nmse_dB;
average_selected_pilots = data.average_selected_pilots;
cfg = data.cfg;

styles = {'^:', 'd-.', 'o-', 'x--', 's--'};
method_colors = { ...
    [0.35, 0.35, 0.35], ...
    [0.4940, 0.1840, 0.5560], ...
    [0.0000, 0.4470, 0.7410], ...
    [0.10, 0.10, 0.10], ...
    [0.8500, 0.3250, 0.0980]};
num_methods = numel(method_names);

set(groot, ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultTextFontName', 'Times New Roman', ...
    'defaultAxesFontSize', 11, ...
    'defaultLineLineWidth', 1.5, ...
    'defaultLineMarkerSize', 6);

file_names = { ...
    'fig_spatial_stability_channel_nmse', ...
    'fig_spatial_stability_pilot_overhead'};

for panel = 1:2
    fig = figure('Color', 'w', 'Position', [100, 100, 450, 340]);
    ax = axes('Parent', fig);
    hold(ax, 'on');
    h = gobjects(num_methods, 1);

    for method = 1:num_methods
        if panel == 1
            values = channel_nmse_dB(:, method);
        else
            values = average_selected_pilots(:, method);
        end
        h(method) = plot(ax, snr_dB, values, styles{method}, ...
            'Color', method_colors{method});
    end

    if panel == 1
        ylabel(ax, 'Channel NMSE (dB)');
        yline(ax, cfg.target_channel_crlb_nmse_dB, 'k:', ...
            'Target', 'HandleVisibility', 'off');
        legend_location = 'southwest';
    else
        ylabel(ax, 'Average pilot symbols');
        ylim(ax, [cfg.min_cycles * cfg.num_reference_states - 2, ...
            cfg.max_cycles * cfg.num_reference_states + 2]);
        legend_location = 'east';
    end

    xlabel(ax, 'Training SNR (dB)');
    grid(ax, 'on');
    box(ax, 'on');
    legend(ax, h, method_names, 'Location', legend_location, ...
        'NumColumns', 1, 'FontSize', 9);

    exportgraphics(fig, fullfile(result_dir, ...
        [file_names{panel}, '.eps']), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(result_dir, ...
        [file_names{panel}, '.png']), 'Resolution', 300);
    savefig(fig, fullfile(result_dir, [file_names{panel}, '.fig']));
    close(fig);
end

fprintf('Replotted Fig. 5 from saved results without simulation.\n');
