%% Replot the saved modulation-order results without rerunning Monte Carlo trials
clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');
data = load(fullfile(result_dir, 'spatial_modulation_order_results.mat'));

colors = [0.0000, 0.4470, 0.7410; ...
          0.8500, 0.3250, 0.0980; ...
          0.4660, 0.6740, 0.1880; ...
          0.4940, 0.1840, 0.5560];
line_styles = {'-', '--', '-.', ':'};
markers = {'o', 's', 'd', '^'};

fig = figure('Color', 'w', 'Position', [100, 100, 520, 400], ...
    'Visible', 'off');
set(groot, 'defaultAxesFontName', 'Times New Roman', ...
    'defaultTextFontName', 'Times New Roman', ...
    'defaultAxesFontSize', 11, 'defaultAxesLineWidth', 0.8, ...
    'defaultLineLineWidth', 1.5, 'defaultLineMarkerSize', 5.5);
hold on;
for i_modulation = 1:numel(data.modulation_orders)
    semilogy(data.snr_dB_vec, data.ser_plot(:, i_modulation), ...
        'Color', colors(i_modulation, :), ...
        'LineStyle', line_styles{i_modulation}, ...
        'Marker', markers{i_modulation}, 'MarkerFaceColor', 'w');
end
set(gca, 'YScale', 'log');
grid on; box on;
xlabel('SNR (dB)');
ylabel('ML SER');
xlim([data.snr_dB_vec(1), data.snr_dB_vec(end)]);
ylim([4e-5, 1]);
xticks(-12:6:24);
legend(data.modulation_names, 'Location', 'southwest', 'FontSize', 9);

set(fig, 'PaperPositionMode', 'auto');
print(fig, fullfile(result_dir, 'fig_spatial_modulation_order_ser.eps'), ...
    '-depsc', '-painters');
print(fig, fullfile(result_dir, 'fig_spatial_modulation_order_ser.png'), ...
    '-dpng', '-r300');
savefig(fig, fullfile(result_dir, 'fig_spatial_modulation_order_ser.fig'));
close(fig);
