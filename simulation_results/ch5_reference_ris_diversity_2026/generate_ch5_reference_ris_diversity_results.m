%% Chapter V-E reference and RIS measurement-diversity simulations
% This self-contained driver consolidates the measurement-diversity tests
% required by Chapter V-E.  Its detector and channel-estimation models are
% adapted from the latest verified Chapter V-C/V-D drivers and, ultimately,
% from tx_and_rx_test4_1.m and tx_and_rx_test5_2.m in the simulation archive.
%
% Observation model:
%   z = |G*s+b|.^2 + w,  w ~ N(0,sigma_w^2 I).
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

rng(20260831, 'twister');

%% Common configuration
cfg.num_users = 2;
cfg.ris_rows = 8;
cfg.ris_cols = 8;
cfg.num_ris = cfg.ris_rows * cfg.ris_cols;
cfg.lambda = 1;
cfg.spacing = cfg.lambda / 2;
cfg.reference_amplitude = 1.5;
cfg.theta_true_deg = [15; 35];
cfg.phi_true_deg = [5; 5];
cfg.path_gains = [exp(1j * 0.35); 0.85 * exp(-1j * 0.55)];
cfg.num_reference_states = 4;

cfg.ris_sweep_snr_dB = -12;
cfg.reference_sweep_snr_dB = -12;
cfg.num_diversity_trials = 10000;
cfg.num_joint_trials = 30000;

cfg.angle_training_snr_dB = 12;
cfg.angle_data_snr_dB = -6;
cfg.angle_num_pilots = 16;
cfg.angle_num_channel_trials = 24;
cfg.angle_data_trials_per_channel = 400;
cfg.angle_perfect_csi_trials = 10000;
cfg.channel_gn_num_restarts = 2;
cfg.channel_gn_max_iter = 30;
cfg.channel_gn_tolerance = 1e-7;
cfg.channel_gn_initial_damping = 1e-2;

qpsk = [1 + 1j; -1 + 1j; -1 - 1j; 1 - 1j] / sqrt(2);
[combination_indices, symbol_candidates] = ...
    enumerate_symbol_vectors(qpsk, cfg.num_users);

fprintf('\nChapter V-E reference and RIS diversity simulations\n');
fprintf('Baseline RIS: %d x %d, users: %d, reference amplitude: %.2f\n', ...
    cfg.ris_rows, cfg.ris_cols, cfg.num_users, cfg.reference_amplitude);

%% 1) Number of RIS intensity observations
ris_side_length = (2:8).';
num_ris_sweep = ris_side_length.^2;
num_ris_cases = numel(ris_side_length);
ris_normalized_dmin = zeros(num_ris_cases, 1);
ris_ml_ser = zeros(num_ris_cases, 1);

for case_index = 1:num_ris_cases
    cfg_case = cfg;
    cfg_case.ris_rows = ris_side_length(case_index);
    cfg_case.ris_cols = ris_side_length(case_index);
    cfg_case.num_ris = num_ris_sweep(case_index);
    H_case = generate_channel_matrix(cfg_case);
    [G_case, b_case] = repeated_reference_measurement( ...
        H_case, cfg.num_reference_states, cfg.reference_amplitude);
    mu_case = intensity_codebook(G_case, b_case, symbol_candidates);
    ris_normalized_dmin(case_index) = normalized_minimum_distance(mu_case);
    ris_ml_ser(case_index) = simulate_ml_ser( ...
        mu_case, combination_indices, cfg.ris_sweep_snr_dB, ...
        cfg.num_diversity_trials, 11000 + case_index);

    fprintf('RIS %2d elements: normalized d_min %.4f, SER %.4e\n', ...
        num_ris_sweep(case_index), ris_normalized_dmin(case_index), ...
        ris_ml_ser(case_index));
end

ris_table = table(num_ris_sweep, ris_normalized_dmin, ris_ml_ser, ...
    'VariableNames', {'Number_of_RIS_elements', 'Normalized_dmin_E', ...
    'ML_SER'});
writetable(ris_table, fullfile(result_dir, 'ris_element_diversity.csv'));

fig = publication_figure([100, 100, 560, 650]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(num_ris_sweep, ris_normalized_dmin, 'o-', ...
    'MarkerFaceColor', [0.25, 0.52, 0.72]);
grid on;
box on;
xlabel('Number of RIS elements, N_r');
ylabel('Normalized d_{min}^{(E)}');
title('(a) Energy-domain separation');
nexttile;
semilogy(num_ris_sweep, plotting_ser(ris_ml_ser, ...
    cfg.num_diversity_trials, cfg.num_users), 's-', ...
    'MarkerFaceColor', [0.85, 0.33, 0.10]);
grid on;
box on;
xlabel('Number of RIS elements, N_r');
ylabel('ML SER');
title(sprintf('(b) Detection at SNR = %g dB', cfg.ris_sweep_snr_dB));
export_publication_figure(fig, result_dir, ...
    'fig_ris_element_measurement_diversity');

%% 2) Number of phase-shifted reference states
H_full = generate_channel_matrix(cfg);
reference_state_sweep = (1:8).';
num_reference_cases = numel(reference_state_sweep);
reference_normalized_dmin = zeros(num_reference_cases, 1);
reference_ml_ser = zeros(num_reference_cases, 1);

for case_index = 1:num_reference_cases
    num_states = reference_state_sweep(case_index);
    [G_case, b_case] = repeated_reference_measurement( ...
        H_full, num_states, cfg.reference_amplitude);
    mu_case = intensity_codebook(G_case, b_case, symbol_candidates);
    reference_normalized_dmin(case_index) = ...
        normalized_minimum_distance(mu_case);
    reference_ml_ser(case_index) = simulate_ml_ser( ...
        mu_case, combination_indices, cfg.reference_sweep_snr_dB, ...
        cfg.num_diversity_trials, 12000 + case_index);

    fprintf('Reference states %d: normalized d_min %.4f, SER %.4e\n', ...
        num_states, reference_normalized_dmin(case_index), ...
        reference_ml_ser(case_index));
end

reference_table = table(reference_state_sweep, ...
    reference_normalized_dmin, reference_ml_ser, ...
    'VariableNames', {'Number_of_reference_states', ...
    'Normalized_dmin_E', 'ML_SER'});
writetable(reference_table, ...
    fullfile(result_dir, 'reference_state_diversity.csv'));

fig = publication_figure([100, 100, 560, 650]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(reference_state_sweep, reference_normalized_dmin, 'o-', ...
    'MarkerFaceColor', [0.25, 0.52, 0.72]);
grid on;
box on;
xlabel('Number of reference states, R');
ylabel('Normalized d_{min}^{(E)}');
title('(a) Energy-domain separation');
nexttile;
semilogy(reference_state_sweep, plotting_ser(reference_ml_ser, ...
    cfg.num_diversity_trials, cfg.num_users), 's-', ...
    'MarkerFaceColor', [0.85, 0.33, 0.10]);
grid on;
box on;
xlabel('Number of reference states, R');
ylabel('ML SER');
title(sprintf('(b) Detection at SNR = %g dB', ...
    cfg.reference_sweep_snr_dB));
export_publication_figure(fig, result_dir, ...
    'fig_reference_state_diversity');

%% 3) User angular separation and end-to-end detection
angle_separation_deg = [0.25; 0.5; 1; 2; 3; 5; 8; 12; 16; 20];
num_angle_cases = numel(angle_separation_deg);
steering_correlation = zeros(num_angle_cases, 1);
normalized_sigma_min = zeros(num_angle_cases, 1);
channel_nmse_dB = zeros(num_angle_cases, 1);
perfect_csi_ser = zeros(num_angle_cases, 1);
estimated_csi_ser = zeros(num_angle_cases, 1);

pilot_symbols = generate_qpsk_pilots( ...
    cfg.angle_num_pilots, cfg.num_users, 13001);
training_reference = cfg.reference_amplitude * ...
    make_uniform_training_reference_phase( ...
    cfg.angle_num_pilots, cfg.num_ris);

for case_index = 1:num_angle_cases
    cfg_angle = cfg;
    cfg_angle.theta_true_deg = [15; 15 + angle_separation_deg(case_index)];
    cfg_angle.phi_true_deg = [5; 5];
    H_angle = generate_channel_matrix(cfg_angle);
    V_angle = H_angle ./ reshape(cfg.path_gains.', 1, []);
    V_normalized = V_angle / sqrt(cfg.num_ris);
    singular_values = svd(V_normalized, 'econ');
    steering_correlation(case_index) = abs( ...
        V_normalized(:, 1)' * V_normalized(:, 2));
    normalized_sigma_min(case_index) = min(singular_values);

    [G_true, b_data] = repeated_reference_measurement( ...
        H_angle, cfg.num_reference_states, cfg.reference_amplitude);
    mu_true = intensity_codebook(G_true, b_data, symbol_candidates);
    perfect_csi_ser(case_index) = simulate_ml_ser( ...
        mu_true, combination_indices, cfg.angle_data_snr_dB, ...
        cfg.angle_perfect_csi_trials, 14000 + case_index);

    nmse_trials = zeros(cfg.angle_num_channel_trials, 1);
    ser_trials = zeros(cfg.angle_num_channel_trials, 1);
    for trial = 1:cfg.angle_num_channel_trials
        H_hat = estimate_channel_realization( ...
            pilot_symbols, training_reference, H_angle, ...
            cfg.angle_training_snr_dB, ...
            15000 + 100 * case_index + trial, cfg);
        nmse_trials(trial) = channel_nmse(H_hat, H_angle);
        G_hat = repmat(H_hat, cfg.num_reference_states, 1);
        ser_trials(trial) = simulate_ml_mismatched_ser( ...
            G_true, G_hat, b_data, combination_indices, ...
            symbol_candidates, cfg.angle_data_snr_dB, ...
            cfg.angle_data_trials_per_channel, ...
            17000 + 100 * case_index + trial);
    end
    channel_nmse_dB(case_index) = 10 * log10(mean(nmse_trials));
    estimated_csi_ser(case_index) = mean(ser_trials);

    fprintf(['Angle separation %5.2f deg: corr %.4f, sigma_min %.4f, ', ...
        'NMSE %.2f dB, perfect/estimated SER %.4e / %.4e\n'], ...
        angle_separation_deg(case_index), steering_correlation(case_index), ...
        normalized_sigma_min(case_index), channel_nmse_dB(case_index), ...
        perfect_csi_ser(case_index), estimated_csi_ser(case_index));
end

angle_table = table(angle_separation_deg, steering_correlation, ...
    normalized_sigma_min, channel_nmse_dB, perfect_csi_ser, ...
    estimated_csi_ser, 'VariableNames', {'Angular_separation_deg', ...
    'Steering_correlation', 'Normalized_sigma_min_V', ...
    'Channel_NMSE_dB', 'Perfect_CSI_ML_SER', 'Estimated_CSI_ML_SER'});
writetable(angle_table, ...
    fullfile(result_dir, 'angular_separation_end_to_end.csv'));

fig = publication_figure([100, 100, 560, 650]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(angle_separation_deg, steering_correlation, 'o-', ...
    'MarkerFaceColor', [0.25, 0.52, 0.72]);
hold on;
plot(angle_separation_deg, normalized_sigma_min, 's--', ...
    'MarkerFaceColor', [0.85, 0.33, 0.10]);
grid on;
box on;
xlabel('User angular separation (deg)');
ylabel('Normalized manifold metric');
legend('Steering-vector correlation', ...
    '\sigma_{min}(V/\surdN_r)', 'Location', 'east');
ylim([0, 1.05]);
title('(a) Spatial conditioning');
nexttile;
semilogy(angle_separation_deg, plotting_ser(perfect_csi_ser, ...
    cfg.angle_perfect_csi_trials, cfg.num_users), 'o-');
hold on;
semilogy(angle_separation_deg, plotting_ser(estimated_csi_ser, ...
    cfg.angle_num_channel_trials * cfg.angle_data_trials_per_channel, ...
    cfg.num_users), 's--');
grid on;
box on;
xlabel('User angular separation (deg)');
ylabel('ML SER');
legend('Perfect CSI', 'Estimated CSI', 'Location', 'northeast');
title(sprintf('(b) End-to-end detection at data SNR = %g dB', ...
    cfg.angle_data_snr_dB));
export_publication_figure(fig, result_dir, ...
    'fig_angular_separation_end_to_end');

%% 4) Reference-only versus joint reference-RIS diversity
% Four states and sixteen element-level observations per state are used in
% both cases.  Reference-only diversity repeats the same central 4x4
% subarray, whereas joint diversity cycles through four interleaved spatial
% groups spanning the full 8x8 aperture.  Thus the total observation count
% and reference power are identical.
reference_only_indices = central_subarray_indices( ...
    cfg.ris_rows, cfg.ris_cols, 4, 4);
joint_index_groups = interleaved_index_groups(cfg.ris_rows, cfg.ris_cols);
state_phases = [0; pi / 2; pi; 3 * pi / 2];
joint_index_groups = select_joint_index_groups( ...
    H_full, symbol_candidates, state_phases, cfg.reference_amplitude, ...
    4096, 18001, joint_index_groups);

G_reference_only = zeros(64, cfg.num_users);
G_joint = zeros(64, cfg.num_users);
b_diversity = zeros(64, 1);
for state = 1:4
    rows = (state - 1) * 16 + (1:16);
    G_reference_only(rows, :) = H_full(reference_only_indices, :);
    G_joint(rows, :) = H_full(joint_index_groups{state}, :);
    b_diversity(rows) = cfg.reference_amplitude * exp(1j * state_phases(state));
end

mu_reference_only = intensity_codebook( ...
    G_reference_only, b_diversity, symbol_candidates);
mu_joint = intensity_codebook(G_joint, b_diversity, symbol_candidates);
joint_design_names = {'Reference only'; 'Joint reference-RIS'};
joint_normalized_dmin = [ ...
    normalized_minimum_distance(mu_reference_only); ...
    normalized_minimum_distance(mu_joint)];

joint_snr_dB = (-16:2:0).';
joint_reference_only_ser = zeros(numel(joint_snr_dB), 1);
joint_reference_ris_ser = zeros(numel(joint_snr_dB), 1);
for snr_index = 1:numel(joint_snr_dB)
    joint_reference_only_ser(snr_index) = simulate_ml_ser( ...
        mu_reference_only, combination_indices, joint_snr_dB(snr_index), ...
        cfg.num_joint_trials, 19000 + snr_index);
    joint_reference_ris_ser(snr_index) = simulate_ml_ser( ...
        mu_joint, combination_indices, joint_snr_dB(snr_index), ...
        cfg.num_joint_trials, 19000 + snr_index);
end

joint_distance_table = table(string(joint_design_names), ...
    joint_normalized_dmin, 'VariableNames', {'Diversity_design', ...
    'Normalized_dmin_E'});
writetable(joint_distance_table, ...
    fullfile(result_dir, 'joint_diversity_distance.csv'));
joint_ser_table = table(joint_snr_dB, joint_reference_only_ser, ...
    joint_reference_ris_ser, 'VariableNames', {'SNR_dB', ...
    'Reference_only_SER', 'Joint_reference_RIS_SER'});
writetable(joint_ser_table, ...
    fullfile(result_dir, 'joint_diversity_ser.csv'));

fprintf('Reference-only normalized d_min: %.4f\n', ...
    joint_normalized_dmin(1));
fprintf('Joint reference-RIS normalized d_min: %.4f\n', ...
    joint_normalized_dmin(2));

fig = publication_figure([100, 100, 560, 650]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
bar(joint_normalized_dmin, 0.62, 'FaceColor', [0.25, 0.52, 0.72]);
grid on;
box on;
set(gca, 'XTick', 1:2, 'XTickLabel', joint_design_names, ...
    'XTickLabelRotation', 8);
ylabel('Normalized d_{min}^{(E)}');
title('(a) Equal-overhead energy separation');
nexttile;
semilogy(joint_snr_dB, plotting_ser(joint_reference_only_ser, ...
    cfg.num_joint_trials, cfg.num_users), 'o-');
hold on;
semilogy(joint_snr_dB, plotting_ser(joint_reference_ris_ser, ...
    cfg.num_joint_trials, cfg.num_users), 's--');
grid on;
box on;
xlabel('SNR (dB)');
ylabel('ML SER');
legend(joint_design_names, 'Location', 'southwest');
title('(b) Detection performance');
export_publication_figure(fig, result_dir, ...
    'fig_joint_reference_ris_diversity');

%% Save complete numerical workspace
save(fullfile(result_dir, 'chapter5_reference_ris_diversity_results.mat'), ...
    'cfg', 'ris_side_length', 'num_ris_sweep', ...
    'ris_normalized_dmin', 'ris_ml_ser', 'reference_state_sweep', ...
    'reference_normalized_dmin', 'reference_ml_ser', ...
    'angle_separation_deg', 'steering_correlation', ...
    'normalized_sigma_min', 'channel_nmse_dB', 'perfect_csi_ser', ...
    'estimated_csi_ser', 'joint_design_names', 'joint_normalized_dmin', ...
    'joint_snr_dB', 'joint_reference_only_ser', ...
    'joint_reference_ris_ser', 'H_full');

fprintf('\nAll outputs saved to:\n%s\n', result_dir);

%% Local functions
function H = generate_channel_matrix(cfg)
    H = zeros(cfg.num_ris, cfg.num_users);
    for user = 1:cfg.num_users
        H(:, user) = cfg.path_gains(user) * steering_vector( ...
            cfg.theta_true_deg(user), cfg.phi_true_deg(user), cfg);
    end
end

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

function [indices, symbols] = enumerate_symbol_vectors(constellation, users)
    constellation = constellation(:);
    order = numel(constellation);
    grids = cell(1, users);
    [grids{:}] = ndgrid(1:order);
    indices = zeros(order^users, users);
    for user = 1:users
        indices(:, user) = grids{user}(:);
    end
    symbols = constellation(indices);
end

function [G, b] = repeated_reference_measurement(H, states, amplitude)
    G = repmat(H, states, 1);
    phase_cycle = [0; pi / 2; pi; 3 * pi / 2];
    phases = phase_cycle(mod((0:states - 1).', 4) + 1);
    phase_matrix = exp(1j * phases) * ones(1, size(H, 1));
    b = amplitude * reshape(phase_matrix.', [], 1);
end

function mu = intensity_codebook(G, b, candidates)
    mu = abs(G * candidates.' + b).^2;
end

function value = normalized_minimum_distance(mu)
    value = minimum_energy_distance(mu) / ...
        max(sqrt(mean(mu(:).^2)), eps);
end

function distance = minimum_energy_distance(mu)
    squared_norms = sum(mu.^2, 1);
    distance_squared = squared_norms.' + squared_norms - ...
        2 * real(mu' * mu);
    distance_squared = max(distance_squared, 0);
    distance_squared(1:size(distance_squared, 1) + 1:end) = inf;
    distance = sqrt(min(distance_squared(:)));
end

function ser = simulate_ml_ser(mu, combination_indices, snr_dB, ...
        num_trials, random_seed)
    rng(random_seed, 'twister');
    num_candidates = size(mu, 2);
    num_users = size(combination_indices, 2);
    true_ids = randi(num_candidates, num_trials, 1);
    noise_variance = mean(mu(:).^2) / 10^(snr_dB / 10);
    observations = mu(:, true_ids) + ...
        sqrt(noise_variance) * randn(size(mu, 1), num_trials);
    metrics = sum(mu.^2, 1).' + sum(observations.^2, 1) - ...
        2 * real(mu' * observations);
    [~, detected_ids] = min(metrics, [], 1);
    detected_ids = detected_ids(:);
    errors = combination_indices(detected_ids, :) ~= ...
        combination_indices(true_ids, :);
    ser = sum(errors(:)) / (num_trials * num_users);
end

function S = generate_qpsk_pilots(num_pilots, num_users, seed)
    old_rng = rng;
    rng(seed, 'twister');
    constellation = [1 + 1j; -1 + 1j; -1 - 1j; 1 - 1j] / sqrt(2);
    indices = randi([1, 4], num_pilots, num_users);
    S = constellation(indices);
    rng(old_rng);
end

function phase_reference = make_uniform_training_reference_phase( ...
        num_pilots, num_ris)
    time_phase = 2 * pi * (0:num_pilots - 1).' / num_pilots;
    space_phase = 2 * pi * (0:num_ris - 1) / num_ris;
    phase_reference = exp(1j * (time_phase + space_phase));
end

function H_hat = estimate_channel_realization(S, B, H_true, snr_dB, ...
        random_seed, cfg)
    rng(random_seed, 'twister');
    mu = abs(S * H_true.' + B).^2;
    noise_variance = mean(mu(:).^2) / 10^(snr_dB / 10);
    z = mu + sqrt(noise_variance) * randn(size(mu));
    H_hat = recover_multiuser_channel_intensity_gn(z, S, B, cfg);
end

function H_hat = recover_multiuser_channel_intensity_gn(z, S, B, cfg)
    H_hat = zeros(size(z, 2), cfg.num_users);
    for element = 1:size(z, 2)
        H_hat(element, :) = recover_single_element_channel_gn( ...
            z(:, element), S, B(:, element), cfg).';
    end
end

function h_best = recover_single_element_channel_gn(z, S, b, cfg)
    users = size(S, 2);
    z = z(:);
    b = b(:);
    weighted = conj(b) .* S;
    linear_matrix = 2 * [real(weighted), -imag(weighted)];
    linear_target = z - abs(b).^2;
    regularization = 1e-3 * trace(linear_matrix.' * linear_matrix) / ...
        max(2 * users, 1);

    if norm(linear_matrix, 'fro') > 1e-12
        x_linear = (linear_matrix.' * linear_matrix + ...
            (regularization + 1e-9) * eye(2 * users)) \ ...
            (linear_matrix.' * linear_target);
        h_linear = x_linear(1:users) + 1j * x_linear(users + 1:end);
    else
        h_linear = zeros(users, 1);
    end

    signal_scale = sqrt(max(mean(max(z, 0)), 1e-3) / max(users, 1));
    h_best = h_linear;
    best_cost = inf;
    for restart = 1:cfg.channel_gn_num_restarts
        if restart == 1 && norm(h_linear) > 1e-10
            h = h_linear;
        else
            h = signal_scale * (randn(users, 1) + ...
                1j * randn(users, 1)) / sqrt(2);
        end
        damping = cfg.channel_gn_initial_damping;
        current_cost = channel_intensity_cost(h, S, b, z);

        for iteration = 1:cfg.channel_gn_max_iter
            q = S * h + b;
            residual = abs(q).^2 - z;
            weighted = conj(q) .* S;
            jacobian = 2 * [real(weighted), -imag(weighted)];
            normal_matrix = jacobian.' * jacobian;
            gradient = jacobian.' * residual;
            diagonal_scaling = diag(max(diag(normal_matrix), 1e-9));
            accepted = false;
            increment = zeros(2 * users, 1);

            for damping_trial = 1:10
                system_matrix = normal_matrix + damping * diagonal_scaling + ...
                    1e-10 * eye(2 * users);
                trial_increment = -system_matrix \ gradient;
                h_trial = h + trial_increment(1:users) + ...
                    1j * trial_increment(users + 1:end);
                trial_cost = channel_intensity_cost(h_trial, S, b, z);
                if trial_cost < current_cost
                    h = h_trial;
                    current_cost = trial_cost;
                    increment = trial_increment;
                    damping = max(damping / 3, 1e-10);
                    accepted = true;
                    break;
                end
                damping = min(damping * 10, 1e10);
            end
            if ~accepted
                break;
            end
            if norm(increment) <= cfg.channel_gn_tolerance * ...
                    (norm([real(h); imag(h)]) + 1)
                break;
            end
        end
        if current_cost < best_cost
            best_cost = current_cost;
            h_best = h;
        end
    end
end

function value = channel_intensity_cost(h, S, b, z)
    value = mean((abs(S * h + b).^2 - z).^2);
end

function value = channel_nmse(H_hat, H_true)
    value = norm(H_hat - H_true, 'fro')^2 / ...
        max(norm(H_true, 'fro')^2, eps);
end

function ser = simulate_ml_mismatched_ser(G_true, G_hat, b, ...
        combination_indices, candidates, snr_dB, num_trials, seed)
    rng(seed, 'twister');
    mu_true = abs(G_true * candidates.' + b).^2;
    mu_hat = abs(G_hat * candidates.' + b).^2;
    num_candidates = size(candidates, 1);
    users = size(combination_indices, 2);
    true_ids = randi(num_candidates, num_trials, 1);
    noise_variance = mean(mu_true(:).^2) / 10^(snr_dB / 10);
    z = mu_true(:, true_ids) + sqrt(noise_variance) * ...
        randn(size(mu_true, 1), num_trials);
    metrics = sum(mu_hat.^2, 1).' + sum(z.^2, 1) - ...
        2 * real(mu_hat' * z);
    [~, detected_ids] = min(metrics, [], 1);
    detected_ids = detected_ids(:);
    errors = combination_indices(detected_ids, :) ~= ...
        combination_indices(true_ids, :);
    ser = sum(errors(:)) / (num_trials * users);
end

function indices = central_subarray_indices(rows, cols, sub_rows, sub_cols)
    row_start = floor((rows - sub_rows) / 2) + 1;
    col_start = floor((cols - sub_cols) / 2) + 1;
    [col_index, row_index] = meshgrid( ...
        col_start:col_start + sub_cols - 1, ...
        row_start:row_start + sub_rows - 1);
    indices = sub2ind([rows, cols], row_index(:), col_index(:));
end

function groups = interleaved_index_groups(rows, cols)
    groups = cell(4, 1);
    group_index = 0;
    for row_parity = 1:2
        for col_parity = 1:2
            group_index = group_index + 1;
            [col_index, row_index] = meshgrid( ...
                col_parity:2:cols, row_parity:2:rows);
            groups{group_index} = sub2ind( ...
                [rows, cols], row_index(:), col_index(:));
        end
    end
end

function best_groups = select_joint_index_groups( ...
        H, candidates, state_phases, amplitude, ...
        num_candidates, seed, initial_groups)
    % Select a state-dependent spatial schedule using the same finite-
    % alphabet distance criterion as the reference codebook. Each candidate
    % uses four state-dependent 16-element masks, so the observation count
    % and per-sample reference power remain fixed.
    old_rng = rng;
    rng(seed, 'twister');
    elements = size(H, 1);
    states = numel(state_phases);
    elements_per_state = elements / states;
    best_groups = initial_groups;
    best_score = joint_schedule_score( ...
        H, candidates, state_phases, amplitude, best_groups);

    for candidate_index = 1:num_candidates
        groups = cell(states, 1);
        for state = 1:states
            groups{state} = randperm(elements, elements_per_state).';
        end
        score = joint_schedule_score( ...
            H, candidates, state_phases, amplitude, groups);
        if score > best_score
            best_score = score;
            best_groups = groups;
        end
    end
    rng(old_rng);
end

function score = joint_schedule_score( ...
        H, candidates, state_phases, amplitude, groups)
    states = numel(groups);
    elements_per_state = numel(groups{1});
    G = zeros(states * elements_per_state, size(H, 2));
    b = zeros(states * elements_per_state, 1);
    for state = 1:states
        rows = (state - 1) * elements_per_state + ...
            (1:elements_per_state);
        G(rows, :) = H(groups{state}, :);
        b(rows) = amplitude * exp(1j * state_phases(state));
    end
    score = normalized_minimum_distance( ...
        intensity_codebook(G, b, candidates));
end

function values = plotting_ser(values, num_trials, num_users)
    floor_value = 0.5 / (num_trials * num_users);
    values = max(values, floor_value);
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
