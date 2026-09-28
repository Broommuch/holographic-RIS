%% Replot spatial-diversity and energy-balance results as separate figures
% This script uses only the saved arrays and does not rerun Monte Carlo trials.
clc;
clear;
close all;

script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');
data = load(fullfile(result_dir, 'spatial_diversity_energy_results.mat'));

set(groot, ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultTextFontName', 'Times New Roman', ...
    'defaultAxesFontSize', 11, ...
    'defaultLineLineWidth', 1.5, ...
    'defaultLineMarkerSize', 6);

trialsN = 10000;
trialsSep = 20000;
trials_code = 1500;
num_codes = 24;
balance_snr = [-6, -3, 0];

fig = new_panel();
semilogx(data.Nvec, data.dminN, 'o-', 'Color', [0, 0.447, 0.741]);
grid on;
box on;
xlabel('Number of receiving units, N_r');
ylabel('Normalized d_{min}^{(E)}');
export_panel(fig, result_dir, 'fig_spatial_unit_distance');

fig = new_panel();
loglog(data.Nvec, max(data.serN, 0.5 / (2 * trialsN)), 's--', ...
    'Color', [0.85, 0.325, 0.098]);
grid on;
box on;
xlabel('Number of receiving units, N_r');
ylabel('ML SER');
export_panel(fig, result_dir, 'fig_spatial_unit_ser');

fig = new_panel();
semilogx(data.sep, data.sigmin, 'o-', 'Color', [0, 0.447, 0.741]);
grid on;
box on;
xlabel('Angular separation (degree)');
ylabel('\sigma_{min}(V)/\surdN_r');
export_panel(fig, result_dir, 'fig_spatial_angular_conditioning');

fig = new_panel();
loglog(data.sep, max(data.serSep, 0.5 / (2 * trialsSep)), 's--', ...
    'Color', [0.85, 0.325, 0.098]);
grid on;
box on;
xlabel('Angular separation (degree)');
ylabel('Perfect-CSI ML SER');
export_panel(fig, result_dir, 'fig_spatial_angular_ser');

fig = new_panel();
errorbar(data.rsr_dB, data.dmean, data.dstd, 'o-', ...
    'Color', [0, 0.447, 0.741], 'CapSize', 4);
grid on;
box on;
xlabel('RSR (dB)');
ylabel('Ensemble-mean d_{min}^{(E)}');
export_panel(fig, result_dir, 'fig_spatial_reference_distance');

fig = new_panel();
styles = {'o-', 's--', 'd-.'};
colors = [0, 0.447, 0.741; 0.85, 0.325, 0.098; 0.929, 0.694, 0.125];
for index = 1:numel(balance_snr)
    semilogy(data.rsr_dB, max(data.serRatio(:, index), ...
        0.5 / (2 * trials_code * num_codes)), styles{index}, ...
        'Color', colors(index, :));
    hold on;
end
grid on;
box on;
xlabel('RSR (dB)');
ylabel('ML SER');
legend(compose('SNR = %d dB', balance_snr), ...
    'Location', 'best', 'FontSize', 9);
export_panel(fig, result_dir, 'fig_spatial_reference_ser');

fprintf('Replotted Figs. 9--11 from saved results without simulation.\n');

function fig = new_panel()
    fig = figure('Color', 'w', 'Position', [100, 100, 450, 340]);
end

function export_panel(fig, result_dir, base_name)
    exportgraphics(fig, fullfile(result_dir, [base_name, '.eps']), ...
        'ContentType', 'vector');
    exportgraphics(fig, fullfile(result_dir, [base_name, '.png']), ...
        'Resolution', 300);
    savefig(fig, fullfile(result_dir, [base_name, '.fig']));
    close(fig);
end
