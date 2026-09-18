%% Regenerate the two summary figures from saved Monte Carlo data
clc;
close all;

script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');
S = load(fullfile(result_dir, 'theory_guided_algorithm_results.mat'));

colors.blue = [0.0000, 0.4470, 0.7410];
colors.orange = [0.8500, 0.3250, 0.0980];
colors.gray = [0.35, 0.35, 0.35];
styles = {'^-', 's--', 'o-'};
method_colors = {colors.gray, colors.orange, colors.blue};
method_names = {'Undamped GN', 'Conventional LM', ...
    'Stability-aware GN'};

fig = publication_figure([80, 80, 920, 350]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
for method = 1:3
    plot(S.separation_deg, S.channel_nmse_dB(:, method), ...
        styles{method}, 'Color', method_colors{method});
    hold on;
end
grid on;
box on;
xlabel('Inter-user angular separation (degree)');
ylabel('Channel NMSE (dB)');
title('(a) Estimation accuracy');
legend(method_names, 'Location', 'southeast');
xticks(S.separation_deg);

nexttile;
for method = 1:3
    plot(S.separation_deg, S.success_probability(:, method), ...
        styles{method}, 'Color', method_colors{method});
    hold on;
end
grid on;
box on;
xlabel('Inter-user angular separation (degree)');
ylabel('Successful recovery probability');
title('(b) Basin reliability');
legend(method_names, 'Location', 'southeast');
xticks(S.separation_deg);
ylim([0, 1.05]);
export_publication_figure(fig, result_dir, ...
    'fig_stability_aware_gn_accuracy_convergence');

fig = publication_figure([80, 80, 920, 350]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
num_decisions = S.e2e_cfg.num_channel_trials * ...
    S.e2e_cfg.num_data_trials * S.e2e_cfg.num_users;
semilogy(S.training_snr_dB, plotting_ser(S.ser_e2e_baseline, ...
    num_decisions), 's--', 'Color', colors.orange);
hold on;
semilogy(S.training_snr_dB, plotting_ser(S.ser_e2e_proposed, ...
    num_decisions), 'o-', 'Color', colors.blue);
semilogy(S.training_snr_dB, plotting_ser(S.ser_e2e_perfect_ml, ...
    num_decisions), 'k-.');
grid on;
box on;
xlabel('Training SNR (dB)');
ylabel('End-to-end SER');
title('(a) Detection reliability');
legend('LM--GN + fixed-budget GS', ...
    'Stability-aware GN + certified GS', 'Perfect-CSI ML', ...
    'Location', 'southwest');
ylim([2e-5, 1]);

nexttile;
plot(S.training_snr_dB, S.channel_nmse_baseline_dB, ...
    's--', 'Color', colors.orange);
hold on;
plot(S.training_snr_dB, S.channel_nmse_proposed_dB, ...
    'o-', 'Color', colors.blue);
grid on;
box on;
xlabel('Training SNR (dB)');
ylabel('Channel NMSE (dB)');
title('(b) Channel-estimation accuracy');
legend('Conventional LM--GN', 'Stability-aware GN', ...
    'Location', 'southwest');
export_publication_figure(fig, result_dir, ...
    'fig_theory_guided_end_to_end');

function values = plotting_ser(values, num_decisions)
    values = max(values, 0.5 / num_decisions);
end

function fig = publication_figure(position)
    fig = figure('Color', 'w', 'Position', position, 'Visible', 'off');
    set(groot, 'defaultAxesFontName', 'Times New Roman');
    set(groot, 'defaultAxesFontSize', 10);
    set(groot, 'defaultAxesLineWidth', 0.8);
    set(groot, 'defaultLineLineWidth', 1.5);
    set(groot, 'defaultLineMarkerSize', 6);
end

function export_publication_figure(fig, result_dir, base_name)
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, fullfile(result_dir, [base_name, '.eps']), ...
        '-depsc', '-painters');
    print(fig, fullfile(result_dir, [base_name, '.png']), ...
        '-dpng', '-r300');
    savefig(fig, fullfile(result_dir, [base_name, '.fig']));
    close(fig);
end
