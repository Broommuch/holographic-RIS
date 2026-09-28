%% Replot perfect-CSI detection as two independent figures
% The saved Monte Carlo results are loaded directly; no trial is rerun.
clc;
clear;
close all;

script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');
data = load(fullfile(result_dir, 'projected_wf_detection_results.mat'));

snr_dB = data.snr_dB;
ser_ml = data.ser_ml;
ser_lin = data.ser_lin;
ser_pwf = data.ser_pwf;
ser_gs = data.ser_gs;
union_bound = data.union_bound;
iterations_pwf = data.iterations_pwf;
iterations_gs = data.iterations_gs;
fixed_budget = data.max_iterations * ones(size(snr_dB));

ser_ml(ser_ml == 0) = NaN;
ser_lin(ser_lin == 0) = NaN;
ser_pwf(ser_pwf == 0) = NaN;
ser_gs(ser_gs == 0) = NaN;

set(groot, ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultTextFontName', 'Times New Roman', ...
    'defaultAxesFontSize', 11, ...
    'defaultLineLineWidth', 1.5, ...
    'defaultLineMarkerSize', 6);

fig = figure('Color', 'w', 'Position', [100, 100, 450, 340]);
semilogy(snr_dB, ser_ml, 'o-', 'Color', [0, 0.447, 0.741]);
hold on;
semilogy(snr_dB, ser_lin, '^-', 'Color', [0.466, 0.674, 0.188]);
semilogy(snr_dB, ser_pwf, 'd-.', 'Color', [0.494, 0.184, 0.556]);
semilogy(snr_dB, ser_gs, 's--', 'Color', [0.85, 0.325, 0.098]);
semilogy(snr_dB, union_bound, 'k-.');
grid on;
box on;
xlabel('SNR (dB)');
ylabel('SER');
ylim([1e-5, 1]);
legend('ML', 'Linearized', 'Projected WF', 'Reference-assisted GS', ...
    'Union bound', 'Location', 'southwest', 'FontSize', 9);
export_panel(fig, result_dir, 'fig_projected_wf_detection_accuracy');

fig = figure('Color', 'w', 'Position', [100, 100, 450, 340]);
plot(snr_dB, fixed_budget, 'k--');
hold on;
plot(snr_dB, iterations_pwf, 'd-.', 'Color', [0.494, 0.184, 0.556]);
plot(snr_dB, iterations_gs, 'o-', 'Color', [0, 0.447, 0.741]);
grid on;
box on;
xlabel('SNR (dB)');
ylabel('Average nonlinear iterations');
legend('Fixed budget', 'Projected WF', 'Reference-assisted GS', ...
    'Location', 'southwest', 'FontSize', 9);
export_panel(fig, result_dir, 'fig_projected_wf_detection_iterations');

fprintf('Replotted Fig. 6 from saved results without simulation.\n');

function export_panel(fig, result_dir, base_name)
    exportgraphics(fig, fullfile(result_dir, [base_name, '.eps']), ...
        'ContentType', 'vector');
    exportgraphics(fig, fullfile(result_dir, [base_name, '.png']), ...
        'Resolution', 300);
    savefig(fig, fullfile(result_dir, [base_name, '.fig']));
    close(fig);
end
