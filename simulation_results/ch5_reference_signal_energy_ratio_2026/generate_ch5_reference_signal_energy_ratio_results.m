%% Reference-to-signal energy-ratio stability experiment for Chapter V-E
% This self-contained script reuses the 8-by-8 RIS geometry, square-law
% observation model, four-phase reference, and post-detection Gaussian-noise
% definition used by the other Chapter V simulations.
%
% To isolate the energy-ratio effect, the total average incident energy per
% element is fixed:
%       P_s + P_b = 1,
% while RSR = 10*log10(P_b/P_s) is swept.  This represents a fixed detector
% power/dynamic-range budget.  The four phase-shifted intensities recover the
% in-phase and quadrature components through opposite-state differencing.
%
% Outputs (EPS, PNG, FIG, CSV, and MAT) are written to ./results.

clc;
clear;
close all;

script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

rng(20260917, 'twister');

%% Configuration inherited from the Chapter V reference-diversity setup
cfg.ris_rows = 8;
cfg.ris_cols = 8;
cfg.num_ris = cfg.ris_rows * cfg.ris_cols;
cfg.lambda = 1;
cfg.spacing = cfg.lambda / 2;
cfg.theta_deg = 15;
cfg.phi_deg = 8;
cfg.path_gain = exp(1j * 0.35);
cfg.num_reference_states = 4;
cfg.total_incident_energy = 1;
cfg.num_monte_carlo = 5000;
cfg.chunk_size = 250;
cfg.phase_success_threshold_deg = 10;

rsr_dB_vec = (-24:3:24).';
snr_dB_vec = [15, 20, 25];
rsr_linear = 10.^(rsr_dB_vec / 10);

signal_energy = cfg.total_incident_energy ./ (1 + rsr_linear);
reference_energy = cfg.total_incident_energy - signal_energy;
interference_depth = 2 * sqrt(signal_energy .* reference_energy) ./ ...
    cfg.total_incident_energy;
normalized_phase_information = interference_depth.^2;

unit_field = cfg.path_gain * steering_vector( ...
    cfg.theta_deg, cfg.phi_deg, cfg);
unit_field = unit_field / sqrt(mean(abs(unit_field).^2));
reference_phases = (0:cfg.num_reference_states - 1).' * ...
    2 * pi / cfg.num_reference_states;

num_ratios = numel(rsr_dB_vec);
num_snr = numel(snr_dB_vec);
phase_rmse_deg = zeros(num_ratios, num_snr);
phase_crlb_deg = zeros(num_ratios, num_snr);
phase_success_probability = zeros(num_ratios, num_snr);
noise_variance = zeros(num_ratios, num_snr);

fprintf('\nReference-to-signal energy-ratio stability simulation\n');
fprintf('RIS: %d x %d, four reference phases, Monte Carlo: %d\n', ...
    cfg.ris_rows, cfg.ris_cols, cfg.num_monte_carlo);
fprintf('Fixed incident-energy budget: P_s + P_b = %.1f\n', ...
    cfg.total_incident_energy);

%% Energy-ratio sweep
for ratio_index = 1:num_ratios
    signal_amplitude = sqrt(signal_energy(ratio_index));
    reference_amplitude = sqrt(reference_energy(ratio_index));
    true_field = signal_amplitude * unit_field;

    reference_matrix = reference_amplitude * ...
        exp(1j * reference_phases) * ones(1, cfg.num_ris);
    field_matrix = ones(cfg.num_reference_states, 1) * true_field.';
    noiseless_intensity = abs(field_matrix + reference_matrix).^2;

    for snr_index = 1:num_snr
        noise_variance(ratio_index, snr_index) = ...
            mean(noiseless_intensity(:).^2) / ...
            10^(snr_dB_vec(snr_index) / 10);

        % For four uniformly spaced reference phases, the scalar phase FIM
        % at every RIS element is 8*P_s*P_b/sigma_w^2.
        phase_crlb_deg(ratio_index, snr_index) = ...
            sqrt(noise_variance(ratio_index, snr_index) / ...
            (8 * signal_energy(ratio_index) * ...
            reference_energy(ratio_index))) * 180 / pi;

        squared_error_sum = 0;
        success_count = 0;
        sample_count = 0;
        trials_completed = 0;
        rng(300000 + 1000 * ratio_index + snr_index, 'twister');

        while trials_completed < cfg.num_monte_carlo
            current_chunk = min(cfg.chunk_size, ...
                cfg.num_monte_carlo - trials_completed);
            observations = noiseless_intensity + ...
                sqrt(noise_variance(ratio_index, snr_index)) * ...
                randn(cfg.num_reference_states, cfg.num_ris, current_chunk);

            real_estimate = reshape( ...
                (observations(1, :, :) - observations(3, :, :)) / ...
                (4 * reference_amplitude), cfg.num_ris, current_chunk);
            imag_estimate = reshape( ...
                (observations(2, :, :) - observations(4, :, :)) / ...
                (4 * reference_amplitude), cfg.num_ris, current_chunk);
            field_estimate = real_estimate + 1j * imag_estimate;

            phase_error = angle(field_estimate .* conj(true_field));
            squared_error_sum = squared_error_sum + ...
                sum(phase_error(:).^2);
            success_count = success_count + sum( ...
                abs(phase_error(:)) <= ...
                cfg.phase_success_threshold_deg * pi / 180);
            sample_count = sample_count + numel(phase_error);
            trials_completed = trials_completed + current_chunk;
        end

        phase_rmse_deg(ratio_index, snr_index) = ...
            sqrt(squared_error_sum / sample_count) * 180 / pi;
        phase_success_probability(ratio_index, snr_index) = ...
            success_count / sample_count;
    end

    fprintf(['RSR %+5.1f dB: kappa %.3f, phase RMSE ', ...
        '[%.2f, %.2f, %.2f] deg\n'], rsr_dB_vec(ratio_index), ...
        interference_depth(ratio_index), phase_rmse_deg(ratio_index, :));
end

%% Save numerical results
summary_table = table(rsr_dB_vec, rsr_linear, signal_energy, ...
    reference_energy, interference_depth, normalized_phase_information, ...
    'VariableNames', {'RSR_dB', 'RSR_linear', 'SignalEnergy', ...
    'ReferenceEnergy', 'InterferenceDepthKappa', ...
    'NormalizedPhaseInformation'});
for snr_index = 1:num_snr
    snr_tag = sprintf('SNR%d', snr_dB_vec(snr_index));
    summary_table.(['PhaseRMSE_', snr_tag, '_deg']) = ...
        phase_rmse_deg(:, snr_index);
    summary_table.(['PhaseCRLB_', snr_tag, '_deg']) = ...
        phase_crlb_deg(:, snr_index);
    summary_table.(['PhaseSuccess_', snr_tag]) = ...
        phase_success_probability(:, snr_index);
    summary_table.(['NoiseVariance_', snr_tag]) = ...
        noise_variance(:, snr_index);
end
writetable(summary_table, fullfile(result_dir, ...
    'reference_signal_energy_ratio_summary.csv'));

save(fullfile(result_dir, ...
    'chapter5_reference_signal_energy_ratio_results.mat'), ...
    'cfg', 'rsr_dB_vec', 'rsr_linear', 'snr_dB_vec', ...
    'signal_energy', 'reference_energy', 'interference_depth', ...
    'normalized_phase_information', 'phase_rmse_deg', ...
    'phase_crlb_deg', 'phase_success_probability', 'noise_variance', ...
    'unit_field', 'reference_phases');

%% Publication figure
fig = publication_figure([100, 100, 900, 360]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(rsr_dB_vec, interference_depth, 'o-', ...
    'Color', [0.0000, 0.4470, 0.7410], ...
    'MarkerFaceColor', 'w', 'LineWidth', 1.6, 'MarkerSize', 5.5);
hold on;
plot(rsr_dB_vec, normalized_phase_information, 's--', ...
    'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor', 'w', 'LineWidth', 1.4, 'MarkerSize', 5.0);
grid on;
box on;
xlabel('Reference-to-signal energy ratio (dB)');
ylabel('Normalized stability metric');
xlim([rsr_dB_vec(1), rsr_dB_vec(end)]);
ylim([0, 1.05]);
xticks(-24:12:24);
legend({'Interference depth, \kappa', ...
    'Relative phase information, \kappa^2'}, ...
    'Location', 'south');
title('(a) Phase-sensitive interference');

nexttile;
colors = [0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980; ...
    0.4660, 0.6740, 0.1880];
markers = {'o', 's', '^'};
mc_handles = gobjects(num_snr, 1);
hold on;
for snr_index = 1:num_snr
    mc_handles(snr_index) = semilogy(rsr_dB_vec, ...
        phase_rmse_deg(:, snr_index), ...
        'Color', colors(snr_index, :), 'LineStyle', '-', ...
        'Marker', markers{snr_index}, 'MarkerFaceColor', 'w', ...
        'LineWidth', 1.6, 'MarkerSize', 5.2);
    semilogy(rsr_dB_vec, phase_crlb_deg(:, snr_index), ...
        'Color', colors(snr_index, :), 'LineStyle', '--', ...
        'LineWidth', 1.1);
end
crlb_handle = semilogy(nan, nan, 'k--', 'LineWidth', 1.1);
yline(cfg.phase_success_threshold_deg, ':', ...
    'Color', [0.35, 0.35, 0.35], 'LineWidth', 1.0);
grid on;
box on;
xlabel('Reference-to-signal energy ratio (dB)');
ylabel('Phase RMSE (degree)');
xlim([rsr_dB_vec(1), rsr_dB_vec(end)]);
ylim([1, 150]);
xticks(-24:12:24);
yticks([1, 2, 5, 10, 20, 50, 100]);
legend([mc_handles; crlb_handle], ...
    {'MC, SNR = 15 dB', 'MC, SNR = 20 dB', ...
    'MC, SNR = 25 dB', 'Corresponding CRLB'}, ...
    'Location', 'north');
title('(b) Four-phase field recovery');

export_publication_figure(fig, result_dir, ...
    'fig_reference_signal_energy_ratio_stability');

%% Report the empirical 10-degree recovery intervals
fprintf('\nEmpirical phase-RMSE intervals below %.1f degrees:\n', ...
    cfg.phase_success_threshold_deg);
for snr_index = 1:num_snr
    valid = phase_rmse_deg(:, snr_index) <= ...
        cfg.phase_success_threshold_deg;
    if any(valid)
        fprintf('SNR %2d dB: RSR from %+g to %+g dB\n', ...
            snr_dB_vec(snr_index), min(rsr_dB_vec(valid)), ...
            max(rsr_dB_vec(valid)));
    else
        fprintf('SNR %2d dB: no sampled RSR satisfies the criterion\n', ...
            snr_dB_vec(snr_index));
    end
end
fprintf('All outputs saved to:\n%s\n', result_dir);

%% Local functions
function v = steering_vector(theta_deg, phi_deg, cfg)
    [y_index, z_index] = meshgrid( ...
        0:cfg.ris_cols - 1, 0:cfg.ris_rows - 1);
    y_position = (y_index(:) - (cfg.ris_cols - 1) / 2) * cfg.spacing;
    z_position = (z_index(:) - (cfg.ris_rows - 1) / 2) * cfg.spacing;
    theta = theta_deg * pi / 180;
    phi = phi_deg * pi / 180;
    wavenumber = 2 * pi / cfg.lambda;
    phase = wavenumber * sin(theta) .* ...
        (cos(phi) * y_position + sin(phi) * z_position);
    v = exp(1j * phase);
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
