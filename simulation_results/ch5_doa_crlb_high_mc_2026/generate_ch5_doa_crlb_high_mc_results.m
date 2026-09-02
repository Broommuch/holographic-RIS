%% Chapter V-B high-Monte-Carlo DOA RMSE and CRLB simulation
% This focused script is derived from the verified Chapter V-B driver:
%   ../ch5_channel_estimation_2026/generate_ch5_channel_estimation_results.m
%
% The signal generation and Gauss--Newton recovery are adapted to the
% post-detection real Gaussian intensity-noise model used in main.tex:
%   z = |S H^T + B|^2 + w,  w ~ N(0,sigma_w^2 I).
%
% Only the DOA-RMSE/CRLB SNR sweep is executed. Increasing the Monte Carlo
% count from 30 to 1000 reduces the sampling fluctuations visible in the
% original figure. Outputs are written to ./results in EPS, PNG, FIG, CSV,
% and MAT formats.

clc;
clear;
close all;

script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

rng(20260828, 'twister');

%% Baseline configuration
cfg.num_users = 2;
cfg.ris_rows = 8;
cfg.ris_cols = 8;
cfg.num_ris = cfg.ris_rows * cfg.ris_cols;
cfg.lambda = 1;
cfg.spacing = cfg.lambda / 2;
cfg.num_pilots = 32;
cfg.reference_amplitude = 1.5;
cfg.num_monte_carlo = 1000;
cfg.fixed_snr_dB = 12;

cfg.theta_true_deg = [15; 35];
cfg.phi_true_deg = [8; -10];
cfg.path_gains = [exp(1j * 0.35); 0.85 * exp(-1j * 0.55)];

cfg.gn_max_iter = 45;
cfg.gn_tolerance = 1e-7;
cfg.gn_initial_damping = 1e-2;
cfg.gn_num_restarts = 2;

cfg.theta_grid_deg = 0:1:60;
cfg.phi_grid_deg = -30:1:30;
cfg.local_refine_max_iter = 60;
cfg.reference_codebook_size = 32;

snr_dB_vec = 0:3:24;
pilot_length_vec = [4, 6, 8, 12, 16, 24, 32, 48];
reference_amplitude_vec = [0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3];
reference_names = {'No reference', 'Constant', 'Random', 'Uniform', 'Optimized'};

max_pilots = max([cfg.num_pilots, pilot_length_vec]);
pilot_symbols_full = generate_qpsk_pilots( ...
    max_pilots, cfg.num_users, 20260829);
pilot_symbols = pilot_symbols_full(1:cfg.num_pilots, :);

H_true = generate_channel_matrix(cfg);

[steering_dictionary, dictionary_theta_deg, dictionary_phi_deg] = ...
    build_steering_dictionary(cfg);

reference_phase_optimized = select_optimized_reference_phase( ...
    pilot_symbols, H_true, cfg.reference_codebook_size, 20260830);
reference_baseline = cfg.reference_amplitude * reference_phase_optimized;

fprintf('\nChapter V-B high-Monte-Carlo DOA/CRLB simulation\n');
fprintf('RIS: %d x %d, users: %d, pilots: %d, Monte Carlo: %d\n', ...
    cfg.ris_rows, cfg.ris_cols, cfg.num_users, ...
    cfg.num_pilots, cfg.num_monte_carlo);

%% 1) SNR sweep: channel NMSE, DOA RMSE, and CRLB
num_snr = numel(snr_dB_vec);
channel_nmse_snr_dB = zeros(num_snr, 1);
doa_rmse_snr_deg = zeros(num_snr, 1);
doa_crlb_snr_deg = zeros(num_snr, 1);
doa_rmse_ci_lower_deg = zeros(num_snr, 1);
doa_rmse_ci_upper_deg = zeros(num_snr, 1);
doa_bias_rms_deg = zeros(num_snr, 1);

mu_baseline = abs(pilot_symbols * H_true.' + reference_baseline).^2;

for i_snr = 1:num_snr
    snr_dB = snr_dB_vec(i_snr);
    noise_variance = mean(mu_baseline(:).^2) / 10^(snr_dB / 10);

    nmse_trials = zeros(cfg.num_monte_carlo, 1);
    doa_sq_error = zeros(cfg.num_monte_carlo, cfg.num_users);
    theta_error = zeros(cfg.num_monte_carlo, cfg.num_users);
    phi_error = zeros(cfg.num_monte_carlo, cfg.num_users);

    for mc = 1:cfg.num_monte_carlo
        rng(100000 + 1000 * i_snr + mc, 'twister');
        z_obs = mu_baseline + sqrt(noise_variance) * randn(size(mu_baseline));

        H_hat = recover_multiuser_channel_intensity_gn( ...
            z_obs, pilot_symbols, reference_baseline, cfg);

        nmse_trials(mc) = norm(H_hat - H_true, 'fro')^2 / ...
            max(norm(H_true, 'fro')^2, eps);

        [theta_hat_deg, phi_hat_deg] = fit_channel_to_doa( ...
            H_hat, steering_dictionary, ...
            dictionary_theta_deg, dictionary_phi_deg, cfg);

        theta_error(mc, :) = (theta_hat_deg - cfg.theta_true_deg).';
        phi_error(mc, :) = (phi_hat_deg - cfg.phi_true_deg).';
        doa_sq_error(mc, :) = theta_error(mc, :).^2 + ...
            phi_error(mc, :).^2;
    end

    channel_nmse_snr_dB(i_snr) = ...
        10 * log10(max(mean(nmse_trials), eps));
    doa_rmse_snr_deg(i_snr) = sqrt(mean(doa_sq_error(:)));
    mse_standard_error = std(doa_sq_error(:), 0) / ...
        sqrt(numel(doa_sq_error));
    rmse_standard_error = mse_standard_error / ...
        max(2 * doa_rmse_snr_deg(i_snr), eps);
    doa_rmse_ci_lower_deg(i_snr) = max( ...
        doa_rmse_snr_deg(i_snr) - 1.96 * rmse_standard_error, 0);
    doa_rmse_ci_upper_deg(i_snr) = ...
        doa_rmse_snr_deg(i_snr) + 1.96 * rmse_standard_error;
    theta_bias = mean(theta_error, 1);
    phi_bias = mean(phi_error, 1);
    doa_bias_rms_deg(i_snr) = sqrt(mean(theta_bias.^2 + phi_bias.^2));
    doa_crlb_snr_deg(i_snr) = compute_doa_crlb_rmse( ...
        pilot_symbols, reference_baseline, cfg, noise_variance);

    fprintf('SNR %5.1f dB: channel NMSE %8.3f dB, DOA RMSE %7.3f deg\n', ...
        snr_dB, channel_nmse_snr_dB(i_snr), doa_rmse_snr_deg(i_snr));
end

snr_table = table( ...
    snr_dB_vec(:), channel_nmse_snr_dB, doa_rmse_snr_deg, ...
    doa_rmse_ci_lower_deg, doa_rmse_ci_upper_deg, ...
    doa_bias_rms_deg, doa_crlb_snr_deg, ...
    'VariableNames', {'SNR_dB', 'Channel_NMSE_dB', ...
    'DOA_RMSE_deg', 'DOA_RMSE_CI95_lower_deg', ...
    'DOA_RMSE_CI95_upper_deg', 'DOA_bias_RMS_deg', ...
    'DOA_CRLB_deg'});
writetable(snr_table, fullfile(result_dir, 'snr_sweep.csv'));

fig = publication_figure([100, 100, 540, 390]);
rmse_plot = max(doa_rmse_snr_deg, 1e-3);
lower_error = rmse_plot - max(doa_rmse_ci_lower_deg, 1e-3);
upper_error = max(doa_rmse_ci_upper_deg, 1e-3) - rmse_plot;
errorbar(snr_dB_vec, rmse_plot, lower_error, upper_error, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 6, 'CapSize', 5);
hold on;
semilogy(snr_dB_vec, max(doa_crlb_snr_deg, 1e-3), 's--', ...
    'LineWidth', 1.5, 'MarkerSize', 6);
set(gca, 'YScale', 'log');
grid on;
box on;
xlabel('SNR (dB)');
ylabel('Joint DOA RMSE (degree)');
legend('Gauss--Newton estimate (95% CI)', 'CRLB', ...
    'Location', 'southwest');
xlim([snr_dB_vec(1), snr_dB_vec(end)]);
export_publication_figure(fig, result_dir, 'fig_doa_rmse_crlb_snr_high_mc');

save(fullfile(result_dir, 'chapter5_doa_crlb_high_mc_results.mat'), ...
    'cfg', 'snr_dB_vec', 'channel_nmse_snr_dB', ...
    'doa_rmse_snr_deg', 'doa_rmse_ci_lower_deg', ...
    'doa_rmse_ci_upper_deg', 'doa_bias_rms_deg', ...
    'doa_crlb_snr_deg', 'H_true', 'pilot_symbols', ...
    'reference_baseline');
fprintf('\nHigh-Monte-Carlo outputs saved to:\n%s\n', result_dir);
return;

%% 2) Pilot-length sweep
num_lengths = numel(pilot_length_vec);
channel_nmse_pilot_dB = zeros(num_lengths, 1);
jacobian_sigma_min_pilot = zeros(num_lengths, 1);

for i_length = 1:num_lengths
    num_pilots = pilot_length_vec(i_length);
    S = pilot_symbols_full(1:num_pilots, :);
    phase_reference = select_optimized_reference_phase( ...
        S, H_true, max(12, floor(cfg.reference_codebook_size / 2)), ...
        200000 + num_pilots);
    B = cfg.reference_amplitude * phase_reference;
    mu = abs(S * H_true.' + B).^2;
    noise_variance = mean(mu(:).^2) / 10^(cfg.fixed_snr_dB / 10);

    nmse_trials = zeros(cfg.num_monte_carlo, 1);
    for mc = 1:cfg.num_monte_carlo
        rng(210000 + 1000 * i_length + mc, 'twister');
        z_obs = mu + sqrt(noise_variance) * randn(size(mu));
        H_hat = recover_multiuser_channel_intensity_gn(z_obs, S, B, cfg);
        nmse_trials(mc) = norm(H_hat - H_true, 'fro')^2 / ...
            max(norm(H_true, 'fro')^2, eps);
    end

    channel_nmse_pilot_dB(i_length) = ...
        10 * log10(max(mean(nmse_trials), eps));
    [jacobian_sigma_min_pilot(i_length), ~] = ...
        channel_jacobian_metrics(S, H_true, B);

    fprintf('Pilots %3d: channel NMSE %8.3f dB\n', ...
        num_pilots, channel_nmse_pilot_dB(i_length));
end

pilot_table = table( ...
    pilot_length_vec(:), channel_nmse_pilot_dB, ...
    jacobian_sigma_min_pilot, ...
    'VariableNames', {'PilotLength', 'Channel_NMSE_dB', ...
    'JacobianSigmaMin'});
writetable(pilot_table, fullfile(result_dir, 'pilot_length_sweep.csv'));

fig = publication_figure([100, 100, 540, 390]);
plot(pilot_length_vec, channel_nmse_pilot_dB, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
box on;
xlabel('Pilot length, T_p');
ylabel('Channel NMSE (dB)');
xticks(pilot_length_vec);
export_publication_figure(fig, result_dir, 'fig_channel_nmse_pilot_length');

%% 3) Reference-amplitude sweep with fixed detector-noise variance
num_amplitudes = numel(reference_amplitude_vec);
channel_nmse_reference_power_dB = zeros(num_amplitudes, 1);
jacobian_sigma_min_reference_power = zeros(num_amplitudes, 1);

baseline_noise_variance = mean(mu_baseline(:).^2) / ...
    10^(cfg.fixed_snr_dB / 10);

for i_amp = 1:num_amplitudes
    reference_amplitude = reference_amplitude_vec(i_amp);
    B = reference_amplitude * reference_phase_optimized;
    mu = abs(pilot_symbols * H_true.' + B).^2;

    nmse_trials = zeros(cfg.num_monte_carlo, 1);
    for mc = 1:cfg.num_monte_carlo
        rng(300000 + 1000 * i_amp + mc, 'twister');
        z_obs = mu + sqrt(baseline_noise_variance) * randn(size(mu));
        H_hat = recover_multiuser_channel_intensity_gn( ...
            z_obs, pilot_symbols, B, cfg);
        nmse_trials(mc) = norm(H_hat - H_true, 'fro')^2 / ...
            max(norm(H_true, 'fro')^2, eps);
    end

    channel_nmse_reference_power_dB(i_amp) = ...
        10 * log10(max(mean(nmse_trials), eps));
    [jacobian_sigma_min_reference_power(i_amp), ~] = ...
        channel_jacobian_metrics(pilot_symbols, H_true, B);

    fprintf('Reference amplitude %4.2f: channel NMSE %8.3f dB\n', ...
        reference_amplitude, channel_nmse_reference_power_dB(i_amp));
end

reference_power_table = table( ...
    reference_amplitude_vec(:), channel_nmse_reference_power_dB, ...
    jacobian_sigma_min_reference_power, ...
    'VariableNames', {'ReferenceAmplitude', 'Channel_NMSE_dB', ...
    'JacobianSigmaMin'});
writetable(reference_power_table, ...
    fullfile(result_dir, 'reference_amplitude_sweep.csv'));

fig = publication_figure([100, 100, 540, 390]);
plot(reference_amplitude_vec, channel_nmse_reference_power_dB, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
box on;
xlabel('Reference amplitude, \rho');
ylabel('Channel NMSE (dB)');
export_publication_figure(fig, result_dir, 'fig_channel_nmse_reference_amplitude');

%% 4) Reference-design comparison at fixed SNR and equal power
num_designs = numel(reference_names);
channel_nmse_reference_design_dB = zeros(num_designs, 1);
jacobian_sigma_min_reference_design = zeros(num_designs, 1);
jacobian_condition_reference_design = zeros(num_designs, 1);

reference_cells = cell(num_designs, 1);
reference_cells{1} = zeros(cfg.num_pilots, cfg.num_ris);
reference_cells{2} = cfg.reference_amplitude * ones(cfg.num_pilots, cfg.num_ris);
reference_cells{3} = cfg.reference_amplitude * make_random_reference_phase( ...
    cfg.num_pilots, cfg.num_ris, 20260831);
reference_cells{4} = cfg.reference_amplitude * make_uniform_reference_phase( ...
    cfg.num_pilots, cfg.num_ris);
reference_cells{5} = reference_baseline;

for i_design = 1:num_designs
    B = reference_cells{i_design};
    mu = abs(pilot_symbols * H_true.' + B).^2;
    noise_variance = mean(mu(:).^2) / 10^(cfg.fixed_snr_dB / 10);

    nmse_trials = zeros(cfg.num_monte_carlo, 1);
    for mc = 1:cfg.num_monte_carlo
        rng(400000 + 1000 * i_design + mc, 'twister');
        z_obs = mu + sqrt(noise_variance) * randn(size(mu));
        H_hat = recover_multiuser_channel_intensity_gn( ...
            z_obs, pilot_symbols, B, cfg);
        nmse_trials(mc) = norm(H_hat - H_true, 'fro')^2 / ...
            max(norm(H_true, 'fro')^2, eps);
    end

    channel_nmse_reference_design_dB(i_design) = ...
        10 * log10(max(mean(nmse_trials), eps));
    [jacobian_sigma_min_reference_design(i_design), ...
        jacobian_condition_reference_design(i_design)] = ...
        channel_jacobian_metrics(pilot_symbols, H_true, B);

    fprintf('%-12s: channel NMSE %8.3f dB, sigma_min %.3e\n', ...
        reference_names{i_design}, ...
        channel_nmse_reference_design_dB(i_design), ...
        jacobian_sigma_min_reference_design(i_design));
end

reference_design_table = table( ...
    string(reference_names(:)), channel_nmse_reference_design_dB, ...
    jacobian_sigma_min_reference_design, ...
    jacobian_condition_reference_design, ...
    'VariableNames', {'ReferenceDesign', 'Channel_NMSE_dB', ...
    'JacobianSigmaMin', 'JacobianConditionNumber'});
writetable(reference_design_table, ...
    fullfile(result_dir, 'reference_design_comparison.csv'));

fig = publication_figure([100, 100, 570, 600]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
bar(1:num_designs, channel_nmse_reference_design_dB, 0.65);
grid on;
box on;
ylabel('Channel NMSE (dB)');
set(gca, 'XTick', 1:num_designs, 'XTickLabel', reference_names, ...
    'XTickLabelRotation', 18);
nexttile;
semilogy(1:num_designs, ...
    max(jacobian_sigma_min_reference_design, 1e-12), ...
    'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
box on;
ylabel('Minimum singular value');
set(gca, 'XTick', 1:num_designs, 'XTickLabel', reference_names, ...
    'XTickLabelRotation', 18);
export_publication_figure(fig, result_dir, ...
    'fig_reference_design_nmse_jacobian');

%% Save all numerical results
save(fullfile(result_dir, 'chapter5_channel_estimation_results.mat'), ...
    'cfg', 'snr_dB_vec', 'channel_nmse_snr_dB', ...
    'doa_rmse_snr_deg', 'doa_crlb_snr_deg', ...
    'pilot_length_vec', 'channel_nmse_pilot_dB', ...
    'jacobian_sigma_min_pilot', ...
    'reference_amplitude_vec', 'channel_nmse_reference_power_dB', ...
    'jacobian_sigma_min_reference_power', ...
    'reference_names', 'channel_nmse_reference_design_dB', ...
    'jacobian_sigma_min_reference_design', ...
    'jacobian_condition_reference_design', ...
    'pilot_symbols', 'H_true', 'reference_baseline');

fprintf('\nAll outputs saved to:\n%s\n', result_dir);

%% Local functions
function S = generate_qpsk_pilots(num_pilots, num_users, seed)
    old_rng = rng;
    rng(seed, 'twister');
    constellation = [1 + 1j; -1 + 1j; -1 - 1j; 1 - 1j] / sqrt(2);
    indices = randi([1, 4], num_pilots, num_users);
    S = constellation(indices);
    rng(old_rng);
end

function H = generate_channel_matrix(cfg)
    H = zeros(cfg.num_ris, cfg.num_users);
    for u = 1:cfg.num_users
        v = steering_vector( ...
            cfg.theta_true_deg(u), cfg.phi_true_deg(u), cfg);
        H(:, u) = cfg.path_gains(u) * v;
    end
end

function v = steering_vector(theta_deg, phi_deg, cfg)
    theta = theta_deg * pi / 180;
    phi = phi_deg * pi / 180;
    [y_index, z_index] = meshgrid( ...
        0:cfg.ris_cols - 1, 0:cfg.ris_rows - 1);
    y_position = (y_index(:) - (cfg.ris_cols - 1) / 2) * cfg.spacing;
    z_position = (z_index(:) - (cfg.ris_rows - 1) / 2) * cfg.spacing;
    k0 = 2 * pi / cfg.lambda;
    ky = k0 * sin(theta) * cos(phi);
    kz = k0 * sin(theta) * sin(phi);
    v = exp(1j * (ky * y_position + kz * z_position));
end

function [v, dv_theta, dv_phi] = steering_vector_with_derivatives( ...
        theta_deg, phi_deg, cfg)
    theta = theta_deg * pi / 180;
    phi = phi_deg * pi / 180;
    [y_index, z_index] = meshgrid( ...
        0:cfg.ris_cols - 1, 0:cfg.ris_rows - 1);
    y_position = (y_index(:) - (cfg.ris_cols - 1) / 2) * cfg.spacing;
    z_position = (z_index(:) - (cfg.ris_rows - 1) / 2) * cfg.spacing;
    k0 = 2 * pi / cfg.lambda;
    phase = k0 * sin(theta) .* ...
        (y_position * cos(phi) + z_position * sin(phi));
    v = exp(1j * phase);
    dphase_theta = k0 * cos(theta) .* ...
        (y_position * cos(phi) + z_position * sin(phi));
    dphase_phi = k0 * sin(theta) .* ...
        (-y_position * sin(phi) + z_position * cos(phi));
    dv_theta = 1j * dphase_theta .* v;
    dv_phi = 1j * dphase_phi .* v;
end

function phase_reference = make_random_reference_phase(num_pilots, num_ris, seed)
    old_rng = rng;
    rng(seed, 'twister');
    phase_reference = exp(1j * 2 * pi * rand(num_pilots, num_ris));
    rng(old_rng);
end

function phase_reference = make_uniform_reference_phase(num_pilots, num_ris)
    time_phase = 2 * pi * (0:num_pilots - 1).' / num_pilots;
    space_phase = 2 * pi * (0:num_ris - 1) / num_ris;
    phase_reference = exp(1j * (time_phase + space_phase));
end

function phase_best = select_optimized_reference_phase(S, H, num_candidates, seed)
    num_pilots = size(S, 1);
    num_ris = size(H, 1);
    phase_best = make_uniform_reference_phase(num_pilots, num_ris);
    [best_score, ~] = channel_jacobian_metrics(S, H, phase_best);

    for candidate = 1:num_candidates
        phase_candidate = make_random_reference_phase( ...
            num_pilots, num_ris, seed + candidate);
        [score, ~] = channel_jacobian_metrics(S, H, phase_candidate);
        if score > best_score
            best_score = score;
            phase_best = phase_candidate;
        end
    end
end

function [sigma_min_global, condition_global] = ...
        channel_jacobian_metrics(S, H, B)
    num_ris = size(H, 1);
    sigma_min_global = inf;
    sigma_max_global = 0;

    for m = 1:num_ris
        h = H(m, :).';
        q = S * h + B(:, m);
        weighted_S = bsxfun(@times, conj(q), S);
        J = 2 * [real(weighted_S), -imag(weighted_S)];
        singular_values = svd(J, 'econ');
        sigma_min_global = min(sigma_min_global, min(singular_values));
        sigma_max_global = max(sigma_max_global, max(singular_values));
    end

    condition_global = sigma_max_global / max(sigma_min_global, eps);
end

function H_hat = recover_multiuser_channel_intensity_gn(z_obs, S, B, cfg)
    num_ris = size(z_obs, 2);
    H_hat = zeros(num_ris, cfg.num_users);
    for m = 1:num_ris
        H_hat(m, :) = recover_single_element_intensity_gn( ...
            z_obs(:, m), S, B(:, m), cfg).';
    end
end

function h_best = recover_single_element_intensity_gn(z, S, b, cfg)
    num_users = size(S, 2);
    z = z(:);
    b = b(:);

    weighted_S = bsxfun(@times, conj(b), S);
    B_real = 2 * [real(weighted_S), -imag(weighted_S)];
    linear_target = z - abs(b).^2;
    regularization = 1e-3 * trace(B_real.' * B_real) / max(2 * num_users, 1);

    if norm(B_real, 'fro') > 1e-12
        x_linear = (B_real.' * B_real + ...
            (regularization + 1e-9) * eye(2 * num_users)) \ ...
            (B_real.' * linear_target);
        h_linear = x_linear(1:num_users) + ...
            1j * x_linear(num_users + 1:end);
    else
        h_linear = zeros(num_users, 1);
    end

    signal_scale = sqrt(max(mean(max(z, 0)), 1e-3) / max(num_users, 1));
    h_best = h_linear;
    cost_best = inf;

    for restart = 1:cfg.gn_num_restarts
        if restart == 1 && norm(h_linear) > 1e-10
            h_current = h_linear;
        else
            h_current = signal_scale * ...
                (randn(num_users, 1) + 1j * randn(num_users, 1)) / sqrt(2);
        end

        damping = cfg.gn_initial_damping;
        cost_current = intensity_cost(h_current, S, b, z);

        for iteration = 1:cfg.gn_max_iter
            q = S * h_current + b;
            residual = abs(q).^2 - z;
            weighted = bsxfun(@times, conj(q), S);
            J = 2 * [real(weighted), -imag(weighted)];
            normal_matrix = J.' * J;
            gradient = J.' * residual;
            diagonal_scaling = diag(max(diag(normal_matrix), 1e-9));
            accepted = false;
            delta = zeros(2 * num_users, 1);

            for damping_trial = 1:10
                system_matrix = normal_matrix + ...
                    damping * diagonal_scaling + 1e-10 * eye(2 * num_users);
                delta_candidate = -system_matrix \ gradient;
                h_candidate = h_current + delta_candidate(1:num_users) + ...
                    1j * delta_candidate(num_users + 1:end);
                cost_candidate = intensity_cost(h_candidate, S, b, z);

                if cost_candidate < cost_current
                    h_current = h_candidate;
                    cost_current = cost_candidate;
                    delta = delta_candidate;
                    damping = max(damping / 3, 1e-10);
                    accepted = true;
                    break;
                end
                damping = min(damping * 10, 1e10);
            end

            if ~accepted
                break;
            end

            if norm(delta) <= cfg.gn_tolerance * ...
                    (norm([real(h_current); imag(h_current)]) + 1)
                break;
            end
        end

        if cost_current < cost_best
            cost_best = cost_current;
            h_best = h_current;
        end
    end
end

function value = intensity_cost(h, S, b, z)
    residual = abs(S * h + b).^2 - z;
    value = mean(residual.^2);
end

function [dictionary, theta_list_deg, phi_list_deg] = ...
        build_steering_dictionary(cfg)
    num_grid = numel(cfg.theta_grid_deg) * numel(cfg.phi_grid_deg);
    dictionary = zeros(cfg.num_ris, num_grid);
    theta_list_deg = zeros(num_grid, 1);
    phi_list_deg = zeros(num_grid, 1);
    index = 0;

    for theta_deg = cfg.theta_grid_deg
        for phi_deg = cfg.phi_grid_deg
            index = index + 1;
            v = steering_vector(theta_deg, phi_deg, cfg);
            dictionary(:, index) = v / max(norm(v), eps);
            theta_list_deg(index) = theta_deg;
            phi_list_deg(index) = phi_deg;
        end
    end
end

function [theta_hat_deg, phi_hat_deg] = fit_channel_to_doa( ...
        H_hat, dictionary, theta_list_deg, phi_list_deg, cfg)
    theta_hat_deg = zeros(cfg.num_users, 1);
    phi_hat_deg = zeros(cfg.num_users, 1);
    options = optimset('Display', 'off', ...
        'MaxIter', cfg.local_refine_max_iter, ...
        'MaxFunEvals', 300, 'TolX', 1e-6, 'TolFun', 1e-10);

    for u = 1:cfg.num_users
        h = H_hat(:, u);
        scores = abs(dictionary' * h).^2 / max(norm(h)^2, eps);
        [~, best_index] = max(scores);
        initial = [theta_list_deg(best_index); phi_list_deg(best_index)];
        objective = @(x) doa_projection_cost(x, h, cfg);
        estimate = fminsearch(objective, initial, options);
        theta_hat_deg(u) = min(max(estimate(1), ...
            min(cfg.theta_grid_deg)), max(cfg.theta_grid_deg));
        phi_hat_deg(u) = min(max(estimate(2), ...
            min(cfg.phi_grid_deg)), max(cfg.phi_grid_deg));
    end
end

function cost = doa_projection_cost(angle_deg, h, cfg)
    theta_min = min(cfg.theta_grid_deg);
    theta_max = max(cfg.theta_grid_deg);
    phi_min = min(cfg.phi_grid_deg);
    phi_max = max(cfg.phi_grid_deg);
    penalty = 1e3 * (max(theta_min - angle_deg(1), 0)^2 + ...
        max(angle_deg(1) - theta_max, 0)^2 + ...
        max(phi_min - angle_deg(2), 0)^2 + ...
        max(angle_deg(2) - phi_max, 0)^2);
    theta_deg = min(max(angle_deg(1), theta_min), theta_max);
    phi_deg = min(max(angle_deg(2), phi_min), phi_max);
    v = steering_vector(theta_deg, phi_deg, cfg);
    alpha = (v' * h) / max(real(v' * v), eps);
    cost = norm(h - alpha * v)^2 / max(norm(h)^2, eps) + penalty;
end

function crlb_rmse_deg = compute_doa_crlb_rmse(S, B, cfg, noise_variance)
    num_users = cfg.num_users;
    num_angle_parameters = 2 * num_users;
    q = S * generate_channel_matrix(cfg).' + B;
    J = zeros(numel(q), 4 * num_users);

    for u = 1:num_users
        [v, dv_theta, dv_phi] = steering_vector_with_derivatives( ...
            cfg.theta_true_deg(u), cfg.phi_true_deg(u), cfg);
        alpha = cfg.path_gains(u);
        derivatives = {alpha * dv_theta, alpha * dv_phi, v, 1j * v};
        columns = [u, num_users + u, ...
            2 * num_users + u, 3 * num_users + u];

        for k = 1:4
            dg = S(:, u) * derivatives{k}.';
            dmu = 2 * real(conj(q) .* dg);
            J(:, columns(k)) = dmu(:);
        end
    end

    F = (J.' * J) / noise_variance;
    F_aa = F(1:num_angle_parameters, 1:num_angle_parameters);
    F_an = F(1:num_angle_parameters, num_angle_parameters + 1:end);
    F_nn = F(num_angle_parameters + 1:end, ...
        num_angle_parameters + 1:end);
    F_doa = F_aa - F_an * pinv(F_nn) * F_an.';
    covariance_bound = pinv(F_doa);
    variances = max(real(diag(covariance_bound)), 0);
    joint_variance = variances(1:num_users) + ...
        variances(num_users + 1:2 * num_users);
    crlb_rmse_deg = sqrt(mean(joint_variance)) * 180 / pi;
end

function fig = publication_figure(position)
    fig = figure('Color', 'w', 'Position', position, 'Visible', 'off');
    set(groot, 'defaultAxesFontName', 'Times New Roman');
    set(groot, 'defaultAxesFontSize', 10);
    set(groot, 'defaultAxesLineWidth', 0.8);
    set(groot, 'defaultLineLineWidth', 1.5);
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
