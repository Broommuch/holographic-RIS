%% Convert the saved joint-per-user DOA metric to a per-angle metric and replot
clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');
mat_file = fullfile(result_dir, 'chapter5_doa_crlb_high_mc_results.mat');
data = load(mat_file);

if ~isfield(data, 'doa_metric_convention') || ...
        ~strcmp(data.doa_metric_convention, 'per-angle')
    scale = 1 / sqrt(2);
    data.doa_rmse_snr_deg = data.doa_rmse_snr_deg * scale;
    data.doa_rmse_ci_lower_deg = data.doa_rmse_ci_lower_deg * scale;
    data.doa_rmse_ci_upper_deg = data.doa_rmse_ci_upper_deg * scale;
    data.doa_bias_rms_deg = data.doa_bias_rms_deg * scale;
    data.doa_crlb_snr_deg = data.doa_crlb_snr_deg * scale;
    data.doa_metric_convention = 'per-angle';
    save(mat_file, '-struct', 'data');
end

snr_table = table(data.snr_dB_vec(:), data.channel_nmse_snr_dB, ...
    data.doa_rmse_snr_deg, data.doa_rmse_ci_lower_deg, ...
    data.doa_rmse_ci_upper_deg, data.doa_bias_rms_deg, ...
    data.doa_crlb_snr_deg, 'VariableNames', {'SNR_dB', ...
    'Channel_NMSE_dB', 'DOA_RMSE_deg', 'DOA_RMSE_CI95_lower_deg', ...
    'DOA_RMSE_CI95_upper_deg', 'DOA_bias_RMS_deg', 'DOA_CRLB_deg'});
writetable(snr_table, fullfile(result_dir, 'snr_sweep.csv'));

fig = figure('Color', 'w', 'Position', [100, 100, 540, 390], ...
    'Visible', 'off');
set(groot, 'defaultAxesFontName', 'Times New Roman', ...
    'defaultTextFontName', 'Times New Roman', ...
    'defaultAxesFontSize', 10, 'defaultAxesLineWidth', 0.8, ...
    'defaultLineLineWidth', 1.5);
rmse_plot = max(data.doa_rmse_snr_deg, 1e-3);
lower_error = rmse_plot - max(data.doa_rmse_ci_lower_deg, 1e-3);
upper_error = max(data.doa_rmse_ci_upper_deg, 1e-3) - rmse_plot;
errorbar(data.snr_dB_vec, rmse_plot, lower_error, upper_error, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 6, 'CapSize', 5);
hold on;
semilogy(data.snr_dB_vec, max(data.doa_crlb_snr_deg, 1e-3), 's--', ...
    'LineWidth', 1.5, 'MarkerSize', 6);
set(gca, 'YScale', 'log');
grid on; box on;
xlabel('SNR (dB)');
ylabel('DOA RMSE (degree)');
legend('Gauss--Newton estimate (95% CI)', 'CRLB', ...
    'Location', 'southwest');
xlim([data.snr_dB_vec(1), data.snr_dB_vec(end)]);

set(fig, 'PaperPositionMode', 'auto');
print(fig, fullfile(result_dir, 'fig_spatial_doa_rmse_crlb.eps'), ...
    '-depsc', '-painters');
print(fig, fullfile(result_dir, 'fig_spatial_doa_rmse_crlb.png'), ...
    '-dpng', '-r300');
savefig(fig, fullfile(result_dir, 'fig_spatial_doa_rmse_crlb.fig'));
close(fig);
