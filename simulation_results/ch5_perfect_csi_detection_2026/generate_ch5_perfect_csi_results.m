%% Chapter V-C perfect-CSI symbol-detection simulations
% This publication-oriented script is based on the newest complete symbol
% detector in the archive:
%   ../01_核心_全息收发联合仿真/02_v2波形域与符号级_2026/tx_and_rx_test4_1.m
%
% The channel, exhaustive ML detector, and reference-aided initialization are
% retained, while the observation noise is changed to the post-detection real
% Gaussian model used in main.tex:
%   z = |G*s+b|^2 + w,  w ~ N(0,sigma_w^2 I).
%
% The script is self-contained and writes EPS, PNG, FIG, CSV, and MAT files to
% ./results.  It covers all results cited in Chapter V-C.

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
cfg.num_reference_states = 4;
cfg.lambda = 1;
cfg.spacing = cfg.lambda / 2;
cfg.reference_amplitude = 1.5;
cfg.theta_true_deg = [15; 35];
cfg.phi_true_deg = [8; -10];
cfg.path_gains = [exp(1j * 0.35); 0.85 * exp(-1j * 0.55)];
cfg.reference_codebook_size = 32;

cfg.gn_max_iter = 35;
cfg.gn_tolerance = 1e-7;
cfg.gn_initial_damping = 1e-2;
cfg.gs_max_iter = 30;
cfg.gs_regularization = 1e-3;

cfg.num_algorithm_trials = 1200;
cfg.num_reference_trials = 1800;
cfg.num_codebook_trials = 3000;
cfg.num_modulation_trials = 1200;
cfg.codebook_test_snr_dB = -6;

snr_algorithm_dB = -15:2:3;
snr_reference_dB = -15:2:3;
snr_modulation_dB = -12:2:8;

qpsk = [1 + 1j; -1 + 1j; -1 - 1j; 1 - 1j] / sqrt(2);
qam16 = generate_square_qam(16);

H_true = generate_channel_matrix(cfg);
G_data = repmat(H_true, cfg.num_reference_states, 1);
[qpsk_indices, qpsk_candidates] = ...
    enumerate_symbol_vectors(qpsk, cfg.num_users);
[qam16_indices, qam16_candidates] = ...
    enumerate_symbol_vectors(qam16, cfg.num_users);

reference_names = {'No reference', 'Constant', 'Random', ...
    'Uniform', 'Optimized'};
[reference_vectors, reference_codebook, optimized_index] = ...
    construct_reference_designs(G_data, qpsk_candidates, cfg);

fprintf('\nChapter V-C perfect-CSI symbol-detection simulations\n');
fprintf('RIS: %d x %d, users: %d, reference states: %d\n', ...
    cfg.ris_rows, cfg.ris_cols, cfg.num_users, ...
    cfg.num_reference_states);
fprintf('Selected reference codeword: %d of %d\n', ...
    optimized_index, cfg.reference_codebook_size);

%% 1) Detector comparison under perfect CSI
b_optimized = reference_vectors(:, 5);
mu_qpsk_optimized = intensity_codebook( ...
    G_data, b_optimized, qpsk_candidates);

num_snr_algorithm = numel(snr_algorithm_dB);
ser_ml = zeros(num_snr_algorithm, 1);
ser_gn = zeros(num_snr_algorithm, 1);
ser_gs = zeros(num_snr_algorithm, 1);

for i_snr = 1:num_snr_algorithm
    [ser_ml(i_snr), ser_gn(i_snr), ser_gs(i_snr)] = ...
        simulate_detector_ser( ...
            G_data, b_optimized, qpsk, qpsk_indices, ...
            qpsk_candidates, mu_qpsk_optimized, ...
            snr_algorithm_dB(i_snr), cfg.num_algorithm_trials, ...
            10000 + i_snr, cfg);

    fprintf(['Detector sweep %5.1f dB: ML %.4e, ', ...
        'projected GN %.4e, reference-assisted GS %.4e\n'], ...
        snr_algorithm_dB(i_snr), ser_ml(i_snr), ...
        ser_gn(i_snr), ser_gs(i_snr));
end

algorithm_table = table(snr_algorithm_dB(:), ser_ml, ser_gn, ser_gs, ...
    'VariableNames', {'SNR_dB', 'ML_SER', 'Projected_GN_SER', ...
    'Reference_Assisted_GS_SER'});
writetable(algorithm_table, fullfile(result_dir, 'algorithm_ser.csv'));

fig = publication_figure([100, 100, 540, 390]);
semilogy(snr_algorithm_dB, plotting_ser(ser_ml, ...
    cfg.num_algorithm_trials, cfg.num_users), 'o-');
hold on;
semilogy(snr_algorithm_dB, plotting_ser(ser_gn, ...
    cfg.num_algorithm_trials, cfg.num_users), 's--');
semilogy(snr_algorithm_dB, plotting_ser(ser_gs, ...
    cfg.num_algorithm_trials, cfg.num_users), 'd-.');
grid on;
box on;
xlabel('SNR (dB)');
ylabel('SER');
legend('ML', 'Projected GN', 'Reference-assisted GS', ...
    'Location', 'southwest');
xlim([snr_algorithm_dB(1), snr_algorithm_dB(end)]);
ylim([2e-4, 1]);
export_publication_figure(fig, result_dir, ...
    'fig_ser_algorithms_perfect_csi');

%% 2) Reference design: ML SER and energy-domain distance
num_reference_designs = numel(reference_names);
num_snr_reference = numel(snr_reference_dB);
ser_reference = zeros(num_snr_reference, num_reference_designs);
dmin_reference = zeros(num_reference_designs, 1);
dmin_reference_normalized = zeros(num_reference_designs, 1);
reference_rms_intensity = zeros(num_reference_designs, 1);

for i_ref = 1:num_reference_designs
    b_current = reference_vectors(:, i_ref);
    mu_current = intensity_codebook(G_data, b_current, qpsk_candidates);
    dmin_reference(i_ref) = minimum_energy_distance(mu_current);
    reference_rms_intensity(i_ref) = sqrt(mean(mu_current(:).^2));
    dmin_reference_normalized(i_ref) = dmin_reference(i_ref) / ...
        max(reference_rms_intensity(i_ref), eps);

    for i_snr = 1:num_snr_reference
        ser_reference(i_snr, i_ref) = simulate_ml_ser( ...
            mu_current, qpsk_indices, snr_reference_dB(i_snr), ...
            cfg.num_reference_trials, 20000 + 100 * i_ref + i_snr);
    end

    fprintf('Reference %-10s: normalized d_min %.4f\n', ...
        reference_names{i_ref}, dmin_reference_normalized(i_ref));
end

reference_ser_table = table(snr_reference_dB(:), ...
    ser_reference(:, 1), ser_reference(:, 2), ser_reference(:, 3), ...
    ser_reference(:, 4), ser_reference(:, 5), ...
    'VariableNames', {'SNR_dB', 'No_reference_SER', 'Constant_SER', ...
    'Random_SER', 'Uniform_SER', 'Optimized_SER'});
writetable(reference_ser_table, ...
    fullfile(result_dir, 'reference_design_ser.csv'));

reference_distance_table = table(reference_names(:), dmin_reference, ...
    dmin_reference_normalized, reference_rms_intensity, ...
    'VariableNames', {'Reference_design', 'dmin_E', ...
    'Normalized_dmin_E', 'RMS_noiseless_intensity'});
writetable(reference_distance_table, ...
    fullfile(result_dir, 'reference_design_distance.csv'));

fig = publication_figure([100, 100, 540, 650]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
styles = {'o-', 's--', '^-.', 'd:', 'p-'};
for i_ref = 1:num_reference_designs
    semilogy(snr_reference_dB, plotting_ser( ...
        ser_reference(:, i_ref), cfg.num_reference_trials, ...
        cfg.num_users), styles{i_ref});
    hold on;
end
grid on;
box on;
xlabel('SNR (dB)');
ylabel('ML SER');
legend(reference_names, 'Location', 'southwest', 'NumColumns', 2);
xlim([snr_reference_dB(1), snr_reference_dB(end)]);
ylim([2e-4, 1]);
title('(a) Detection performance');

nexttile;
bar(dmin_reference_normalized, 0.68, ...
    'FaceColor', [0.25, 0.52, 0.72]);
grid on;
box on;
set(gca, 'XTick', 1:num_reference_designs, ...
    'XTickLabel', reference_names, 'XTickLabelRotation', 18);
ylabel('Normalized d_{min}^{(E)}');
title('(b) Energy-domain separation');
export_publication_figure(fig, result_dir, ...
    'fig_reference_design_perfect_csi');

%% 3) Relationship between codebook distance and SER
num_codewords = size(reference_codebook, 2);
dmin_codebook = zeros(num_codewords, 1);
dmin_codebook_normalized = zeros(num_codewords, 1);
ser_codebook = zeros(num_codewords, 1);

for i_code = 1:num_codewords
    b_current = reference_codebook(:, i_code);
    mu_current = intensity_codebook(G_data, b_current, qpsk_candidates);
    dmin_codebook(i_code) = minimum_energy_distance(mu_current);
    rms_current = sqrt(mean(mu_current(:).^2));
    dmin_codebook_normalized(i_code) = dmin_codebook(i_code) / ...
        max(rms_current, eps);
    ser_codebook(i_code) = simulate_ml_ser( ...
        mu_current, qpsk_indices, cfg.codebook_test_snr_dB, ...
        cfg.num_codebook_trials, 30000 + i_code);
end

ser_codebook_plot = plotting_ser(ser_codebook, ...
    cfg.num_codebook_trials, cfg.num_users);
correlation_matrix = corrcoef(dmin_codebook_normalized, ...
    log10(ser_codebook_plot));
distance_log_ser_correlation = correlation_matrix(1, 2);

codebook_table = table((1:num_codewords).', dmin_codebook, ...
    dmin_codebook_normalized, ser_codebook, ...
    'VariableNames', {'Codeword_index', 'dmin_E', ...
    'Normalized_dmin_E', 'ML_SER'});
writetable(codebook_table, fullfile(result_dir, 'codebook_dmin_ser.csv'));

fig = publication_figure([100, 100, 540, 390]);
semilogy(dmin_codebook_normalized, ser_codebook_plot, 'o', ...
    'MarkerSize', 6, 'MarkerFaceColor', [0.25, 0.52, 0.72]);
hold on;
fit_coefficients = polyfit(dmin_codebook_normalized, ...
    log10(ser_codebook_plot), 1);
fit_x = linspace(min(dmin_codebook_normalized), ...
    max(dmin_codebook_normalized), 100);
fit_y = 10.^polyval(fit_coefficients, fit_x);
semilogy(fit_x, fit_y, '--', 'Color', [0.85, 0.33, 0.10]);
grid on;
box on;
xlabel('Normalized d_{min}^{(E)}');
ylabel('ML SER');
legend('Reference codewords', 'Log-linear trend', ...
    'Location', 'southwest');
text(0.97, 0.95, sprintf('corr. = %.2f', ...
    distance_log_ser_correlation), 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
export_publication_figure(fig, result_dir, 'fig_ser_vs_dmin');

%% 4) Modulation order under the same optimized reference
mu_qam16_optimized = intensity_codebook( ...
    G_data, b_optimized, qam16_candidates);
num_snr_modulation = numel(snr_modulation_dB);
ser_qpsk_modulation = zeros(num_snr_modulation, 1);
ser_qam16_modulation = zeros(num_snr_modulation, 1);

for i_snr = 1:num_snr_modulation
    ser_qpsk_modulation(i_snr) = simulate_ml_ser( ...
        mu_qpsk_optimized, qpsk_indices, snr_modulation_dB(i_snr), ...
        cfg.num_modulation_trials, 40000 + i_snr);
    ser_qam16_modulation(i_snr) = simulate_ml_ser( ...
        mu_qam16_optimized, qam16_indices, snr_modulation_dB(i_snr), ...
        cfg.num_modulation_trials, 50000 + i_snr);

    fprintf('Modulation sweep %5.1f dB: QPSK %.4e, 16-QAM %.4e\n', ...
        snr_modulation_dB(i_snr), ser_qpsk_modulation(i_snr), ...
        ser_qam16_modulation(i_snr));
end

dmin_qpsk = minimum_energy_distance(mu_qpsk_optimized);
dmin_qam16 = minimum_energy_distance(mu_qam16_optimized);
normalized_dmin_qpsk = dmin_qpsk / ...
    sqrt(mean(mu_qpsk_optimized(:).^2));
normalized_dmin_qam16 = dmin_qam16 / ...
    sqrt(mean(mu_qam16_optimized(:).^2));

modulation_table = table(snr_modulation_dB(:), ...
    ser_qpsk_modulation, ser_qam16_modulation, ...
    'VariableNames', {'SNR_dB', 'QPSK_SER', 'QAM16_SER'});
writetable(modulation_table, fullfile(result_dir, 'modulation_ser.csv'));

modulation_distance_table = table({'QPSK'; '16-QAM'}, ...
    [dmin_qpsk; dmin_qam16], ...
    [normalized_dmin_qpsk; normalized_dmin_qam16], ...
    'VariableNames', {'Modulation', 'dmin_E', 'Normalized_dmin_E'});
writetable(modulation_distance_table, ...
    fullfile(result_dir, 'modulation_distance.csv'));

fig = publication_figure([100, 100, 540, 390]);
semilogy(snr_modulation_dB, plotting_ser(ser_qpsk_modulation, ...
    cfg.num_modulation_trials, cfg.num_users), 'o-');
hold on;
semilogy(snr_modulation_dB, plotting_ser(ser_qam16_modulation, ...
    cfg.num_modulation_trials, cfg.num_users), 's--');
grid on;
box on;
xlabel('SNR (dB)');
ylabel('ML SER');
legend('QPSK', '16-QAM', 'Location', 'southwest');
xlim([snr_modulation_dB(1), snr_modulation_dB(end)]);
ylim([2e-4, 1]);
export_publication_figure(fig, result_dir, 'fig_ser_modulation_order');

%% Save complete workspace data
save(fullfile(result_dir, 'chapter5_perfect_csi_detection_results.mat'), ...
    'cfg', 'snr_algorithm_dB', 'ser_ml', 'ser_gn', 'ser_gs', ...
    'snr_reference_dB', 'reference_names', 'ser_reference', ...
    'dmin_reference', 'dmin_reference_normalized', ...
    'reference_rms_intensity', 'optimized_index', ...
    'dmin_codebook', 'dmin_codebook_normalized', 'ser_codebook', ...
    'distance_log_ser_correlation', 'snr_modulation_dB', ...
    'ser_qpsk_modulation', 'ser_qam16_modulation', ...
    'dmin_qpsk', 'dmin_qam16', 'normalized_dmin_qpsk', ...
    'normalized_dmin_qam16', 'H_true', 'reference_vectors', ...
    'reference_codebook');

fprintf('\nAll outputs saved to:\n%s\n', result_dir);

%% Local functions
function H = generate_channel_matrix(cfg)
    H = zeros(cfg.num_ris, cfg.num_users);
    for u = 1:cfg.num_users
        H(:, u) = cfg.path_gains(u) * steering_vector( ...
            cfg.theta_true_deg(u), cfg.phi_true_deg(u), cfg);
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
    phase = wavenumber * (sin(theta) * cos(phi) * y_position + ...
        sin(theta) * sin(phi) * z_position);
    v = exp(1j * phase);
end

function constellation = generate_square_qam(order)
    side = sqrt(order);
    levels = -(side - 1):2:(side - 1);
    [real_grid, imag_grid] = meshgrid(levels, fliplr(levels));
    constellation = real_grid(:) + 1j * imag_grid(:);
    constellation = constellation / sqrt(mean(abs(constellation).^2));
end

function [indices, symbols] = enumerate_symbol_vectors(constellation, users)
    constellation = constellation(:);
    order = numel(constellation);
    grids = cell(1, users);
    [grids{:}] = ndgrid(1:order);
    indices = zeros(order^users, users);
    for u = 1:users
        indices(:, u) = grids{u}(:);
    end
    symbols = constellation(indices);
end

function [reference_vectors, codebook, optimized_index] = ...
        construct_reference_designs(G, candidates, cfg)
    rows = cfg.num_reference_states;
    elements = cfg.num_ris;
    amplitude = cfg.reference_amplitude;

    b_none = zeros(rows * elements, 1);
    b_constant = amplitude * ones(rows * elements, 1);

    rng(20260830, 'twister');
    random_matrix = amplitude * exp(1j * 2 * pi * rand(rows, elements));
    b_random = reshape(random_matrix.', [], 1);

    state_phase = (0:rows - 1).' * 2 * pi / rows;
    uniform_matrix = amplitude * exp(1j * state_phase) * ones(1, elements);
    b_uniform = reshape(uniform_matrix.', [], 1);

    codebook = zeros(rows * elements, cfg.reference_codebook_size);
    codebook(:, 1) = b_uniform;
    codebook(:, 2) = b_constant;
    codebook(:, 3) = b_random;
    phase_alphabet = (0:rows - 1) * 2 * pi / rows;

    num_observations = rows * elements;
    for i_code = 4:cfg.reference_codebook_size
        % Vary both the phase pattern and its spatial/state support while
        % preserving ||b||_2^2 = Q*amplitude^2 for every codeword.  This
        % produces a meaningful range of energy-domain separations without
        % conflating distance with total reference power.
        support_fraction = 0.08 + 0.92 * ...
            (i_code - 4) / max(cfg.reference_codebook_size - 4, 1);
        num_active = max(rows, round(support_fraction * num_observations));
        active_indices = randperm(num_observations, num_active);
        active_amplitude = amplitude * sqrt(num_observations / num_active);
        active_phases = phase_alphabet(randi(rows, num_active, 1));
        codeword = zeros(num_observations, 1);
        codeword(active_indices) = active_amplitude * ...
            exp(1j * active_phases(:));
        codebook(:, i_code) = codeword;
    end

    dmin_codebook = zeros(cfg.reference_codebook_size, 1);
    for i_code = 1:cfg.reference_codebook_size
        mu = intensity_codebook(G, codebook(:, i_code), candidates);
        dmin_codebook(i_code) = minimum_energy_distance(mu) / ...
            max(sqrt(mean(mu(:).^2)), eps);
    end
    [~, optimized_index] = max(dmin_codebook);
    b_optimized = codebook(:, optimized_index);

    reference_vectors = [b_none, b_constant, b_random, ...
        b_uniform, b_optimized];
end

function mu = intensity_codebook(G, b, candidates)
    mu = abs(G * candidates.' + b).^2;
end

function dmin = minimum_energy_distance(mu)
    squared_norms = sum(mu.^2, 1);
    distance_squared = squared_norms.' + squared_norms - ...
        2 * real(mu' * mu);
    distance_squared = max(distance_squared, 0);
    distance_squared(1:size(distance_squared, 1) + 1:end) = inf;
    dmin = sqrt(min(distance_squared(:)));
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

    mu_norms = sum(mu.^2, 1).';
    observation_norms = sum(observations.^2, 1);
    metrics = mu_norms + observation_norms - 2 * real(mu' * observations);
    [~, detected_ids] = min(metrics, [], 1);
    detected_ids = detected_ids(:);

    errors = combination_indices(detected_ids, :) ~= ...
        combination_indices(true_ids, :);
    ser = sum(errors(:)) / (num_trials * num_users);
end

function [ser_ml, ser_gn, ser_gs] = simulate_detector_ser( ...
        G, b, constellation, combination_indices, candidates, mu, ...
        snr_dB, num_trials, random_seed, cfg)
    rng(random_seed, 'twister');
    num_candidates = size(mu, 2);
    num_users = size(combination_indices, 2);
    true_ids = randi(num_candidates, num_trials, 1);
    noise_variance = mean(mu(:).^2) / 10^(snr_dB / 10);
    observations = mu(:, true_ids) + ...
        sqrt(noise_variance) * randn(size(mu, 1), num_trials);

    mu_norms = sum(mu.^2, 1).';
    observation_norms = sum(observations.^2, 1);
    metrics = mu_norms + observation_norms - 2 * real(mu' * observations);
    [~, detected_ml] = min(metrics, [], 1);
    detected_ml = detected_ml(:);

    detected_gn_indices = zeros(num_trials, num_users);
    detected_gs_indices = zeros(num_trials, num_users);
    back_projection = (G' * G + ...
        cfg.gs_regularization * eye(num_users)) \ G';

    for trial = 1:num_trials
        z = observations(:, trial);
        initial = reference_linearized_initialization(z, G, b);
        detected_gn_indices(trial, :) = projected_gn_detector( ...
            z, G, b, constellation, initial, cfg).';
        detected_gs_indices(trial, :) = reference_assisted_gs_detector( ...
            z, G, b, constellation, initial, back_projection, cfg).';
    end

    true_indices = combination_indices(true_ids, :);
    ml_indices = combination_indices(detected_ml, :);
    ser_ml = sum(ml_indices(:) ~= true_indices(:)) / ...
        (num_trials * num_users);
    ser_gn = sum(detected_gn_indices(:) ~= true_indices(:)) / ...
        (num_trials * num_users);
    ser_gs = sum(detected_gs_indices(:) ~= true_indices(:)) / ...
        (num_trials * num_users);

    %#ok<NASGU> candidates is retained in the interface to document the common
    % finite-alphabet search space used by all three detectors.
end

function initial = reference_linearized_initialization(z, G, b)
    users = size(G, 2);
    weighted_matrix = conj(b) .* G;
    linear_matrix = 2 * [real(weighted_matrix), -imag(weighted_matrix)];
    right_hand_side = z - abs(b).^2;
    regularization = 1e-3 * trace(linear_matrix.' * linear_matrix) / ...
        max(2 * users, 1);
    real_initial = (linear_matrix.' * linear_matrix + ...
        regularization * eye(2 * users)) \ ...
        (linear_matrix.' * right_hand_side);
    initial = real_initial(1:users) + ...
        1j * real_initial(users + 1:end);
end

function detected_indices = projected_gn_detector( ...
        z, G, b, constellation, initial, cfg)
    users = size(G, 2);
    initial_points = [initial, zeros(users, 1)];
    best_cost = inf;
    best_indices = ones(users, 1);

    for restart = 1:size(initial_points, 2)
        s = initial_points(:, restart);
        damping = cfg.gn_initial_damping;
        current_cost = norm(z - abs(G * s + b).^2)^2;

        for iteration = 1:cfg.gn_max_iter
            affine_field = G * s + b;
            predicted = abs(affine_field).^2;
            weighted_matrix = conj(affine_field) .* G;
            jacobian = 2 * [real(weighted_matrix), -imag(weighted_matrix)];
            residual = z - predicted;
            increment = (jacobian.' * jacobian + ...
                damping * eye(2 * users)) \ (jacobian.' * residual);
            s_trial = s + increment(1:users) + ...
                1j * increment(users + 1:end);
            trial_cost = norm(z - abs(G * s_trial + b).^2)^2;

            if trial_cost < current_cost
                s = s_trial;
                current_cost = trial_cost;
                damping = max(damping / 2, 1e-9);
            else
                damping = min(damping * 5, 1e9);
            end

            if norm(increment) / max(norm([real(s); imag(s)]), 1) < ...
                    cfg.gn_tolerance
                break;
            end
        end

        indices = quantize_indices(s, constellation);
        s_projected = constellation(indices);
        projected_cost = norm(z - abs(G * s_projected + b).^2)^2;
        if projected_cost < best_cost
            best_cost = projected_cost;
            best_indices = indices;
        end
    end
    detected_indices = best_indices;
end

function detected_indices = reference_assisted_gs_detector( ...
        z, G, b, constellation, initial, back_projection, cfg)
    users = size(G, 2);
    initial_indices = quantize_indices(initial, constellation);
    initial_points = [constellation(initial_indices), zeros(users, 1)];
    measured_magnitude = sqrt(max(z, 0));
    best_cost = inf;
    best_indices = initial_indices;

    for restart = 1:size(initial_points, 2)
        s = initial_points(:, restart);
        for iteration = 1:cfg.gs_max_iter
            affine_field = G * s + b;
            phase = affine_field ./ max(abs(affine_field), eps);
            projected_field = measured_magnitude .* phase;
            continuous_update = back_projection * (projected_field - b);
            indices = quantize_indices(continuous_update, constellation);
            s_next = constellation(indices);
            if norm(s_next - s) / max(norm(s), 1) < cfg.gn_tolerance
                s = s_next;
                break;
            end
            s = s_next;
        end

        indices = quantize_indices(s, constellation);
        s_projected = constellation(indices);
        cost = norm(z - abs(G * s_projected + b).^2)^2;
        if cost < best_cost
            best_cost = cost;
            best_indices = indices;
        end
    end
    detected_indices = best_indices;
end

function indices = quantize_indices(symbols, constellation)
    symbols = symbols(:);
    constellation = constellation(:).';
    distances = abs(symbols - constellation);
    [~, indices] = min(distances, [], 2);
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
