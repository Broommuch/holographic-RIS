%% Algorithm-1 practical-initialization simulations for Chapter V
% This self-contained script repeats the GN and end-to-end experiments
% without access to the true parameters.  It implements the complete
% initialization in Algorithm 1: four-phase field recovery, unstructured
% channel LS, coarse angular search, gain projection, and parametric GN.
% EPS, PNG, FIG, CSV, and MAT outputs are written to ./results.

clc;
clear;
close all;

script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

rng(20260918, 'twister');
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');
quick_check = strcmp(getenv('HOLO_QUICK_CHECK'), '1');

%% Common configuration
cfg.num_users = 2;
cfg.ris_rows = 8;
cfg.ris_cols = 8;
cfg.num_ris = cfg.ris_rows * cfg.ris_cols;
cfg.lambda = 1;
cfg.spacing = cfg.lambda / 2;
cfg.reference_amplitude = 1.5;
cfg.num_reference_states = 4;
cfg.path_gains = [exp(1j * 0.35); 0.85 * exp(-1j * 0.55)];

cfg.gn_max_iter = 30;
cfg.gn_tolerance = 1e-6;
cfg.gn_initial_damping = 1e-3;
cfg.gn_target_condition = 200;
cfg.gs_max_iter = 25;
cfg.gs_tolerance = 1e-7;
cfg.gs_regularization = 1e-3;
cfg.coarse_theta_grid_deg = 0:0.5:60;
cfg.coarse_phi_grid_deg = -30:1:30;
cfg.pilot_shared_fraction = 0.75;

qpsk = [1 + 1j; -1 + 1j; -1 - 1j; 1 - 1j] / sqrt(2);
[candidate_indices, candidate_symbols] = ...
    enumerate_symbol_vectors(qpsk, cfg.num_users);

colors.blue = [0.0000, 0.4470, 0.7410];
colors.orange = [0.8500, 0.3250, 0.0980];
colors.yellow = [0.9290, 0.6940, 0.1250];
colors.purple = [0.4940, 0.1840, 0.5560];
colors.gray = [0.35, 0.35, 0.35];

fprintf('\nAlgorithm-1 practical-initialization simulations\n');
fprintf('RIS: %d x %d, users: %d\n', ...
    cfg.ris_rows, cfg.ris_cols, cfg.num_users);
[coarse_dictionary, coarse_theta_deg, coarse_phi_deg] = ...
    make_coarse_dictionary(cfg);

%% 1) Stability-aware parametric GN
gn_cfg = cfg;
gn_cfg.num_pilots = 16;
gn_cfg.snr_dB = 12;
gn_cfg.num_monte_carlo = 300;
if quick_check
    gn_cfg.num_monte_carlo = 6;
end
gn_cfg.theta_center_deg = 20;
gn_cfg.phi_true_deg = [5; 5];
gn_cfg.success_threshold_deg = 1;

separation_deg = [0.5, 1, 2, 4, 8, 16];
method_names_gn = {'Undamped GN', 'Conventional LM', ...
    'Stability-aware GN'};
num_gn_methods = numel(method_names_gn);
num_separations = numel(separation_deg);

success_probability = zeros(num_separations, num_gn_methods);
doa_rmse_deg = zeros(num_separations, num_gn_methods);
channel_nmse_dB = zeros(num_separations, num_gn_methods);
average_gn_iterations = zeros(num_separations, num_gn_methods);
average_rejections = zeros(num_separations, num_gn_methods);
initial_doa_rmse_deg = zeros(num_separations, 1);
initial_channel_nmse_dB = zeros(num_separations, 1);
trace_store = cell(num_separations, num_gn_methods);

% T_p counts intensity snapshots.  Four phase states are used per pilot
% vector, so T_p=16 provides four field snapshots.  The pilot matrix is
% correlated but full rank, which is required by Algorithm 1's LS step.
[pilot_symbols, training_reference, pilot_cycles] = ...
    make_four_phase_training(gn_cfg, 20260919);
pilot_rank = rank(pilot_cycles);
pilot_condition = cond(pilot_cycles' * pilot_cycles);
fprintf('Training pilot rank %d, Gram condition number %.3f\n', ...
    pilot_rank, pilot_condition);

for i_sep = 1:num_separations
    theta_true_deg = gn_cfg.theta_center_deg + ...
        0.5 * separation_deg(i_sep) * [-1; 1];
    eta_true = pack_eta(theta_true_deg, gn_cfg.phi_true_deg, ...
        gn_cfg.path_gains);
    [mu_true, ~, H_true] = training_model_and_jacobian( ...
        eta_true, pilot_symbols, training_reference, gn_cfg);
    noise_variance = mean(mu_true.^2) / 10^(gn_cfg.snr_dB / 10);

    sq_angle_error = zeros(gn_cfg.num_monte_carlo, num_gn_methods);
    nmse_trial = zeros(gn_cfg.num_monte_carlo, num_gn_methods);
    iteration_trial = zeros(gn_cfg.num_monte_carlo, num_gn_methods);
    rejection_trial = zeros(gn_cfg.num_monte_carlo, num_gn_methods);
    initial_angle_sq_error = zeros(gn_cfg.num_monte_carlo, 1);
    initial_nmse_trial = zeros(gn_cfg.num_monte_carlo, 1);
    cost_trace_trial = nan(gn_cfg.num_monte_carlo, ...
        gn_cfg.gn_max_iter + 1, num_gn_methods);

    for mc = 1:gn_cfg.num_monte_carlo
        rng(100000 + 1000 * i_sep + mc, 'twister');
        z = mu_true + sqrt(noise_variance) * randn(size(mu_true));
        [eta_initial, H_initial] = algorithm1_initialization( ...
            z, pilot_cycles, gn_cfg, coarse_dictionary, ...
            coarse_theta_deg, coarse_phi_deg);
        [~, initial_angle_sq_error(mc)] = matched_doa_error( ...
            eta_initial, theta_true_deg, gn_cfg.phi_true_deg, gn_cfg);
        initial_nmse_trial(mc) = norm(H_initial - H_true, 'fro')^2 / ...
            max(norm(H_true, 'fro')^2, eps);

        for i_method = 1:num_gn_methods
            [eta_hat, information] = parametric_gn( ...
                z, pilot_symbols, training_reference, eta_initial, ...
                i_method, gn_cfg);
            H_hat = eta_to_channel(eta_hat, gn_cfg);
            [angle_rmse, angle_sq_error] = matched_doa_error( ...
                eta_hat, theta_true_deg, gn_cfg.phi_true_deg, gn_cfg);
            sq_angle_error(mc, i_method) = angle_sq_error;
            nmse_trial(mc, i_method) = norm(H_hat - H_true, 'fro')^2 / ...
                max(norm(H_true, 'fro')^2, eps);
            iteration_trial(mc, i_method) = information.iterations;
            rejection_trial(mc, i_method) = information.rejections;
            cost_trace_trial(mc, :, i_method) = information.cost_trace;
            information.success = angle_rmse < gn_cfg.success_threshold_deg; %#ok<NASGU>
        end
    end

    for i_method = 1:num_gn_methods
        per_trial_rmse = sqrt(sq_angle_error(:, i_method));
        success_probability(i_sep, i_method) = mean( ...
            per_trial_rmse < gn_cfg.success_threshold_deg);
        doa_rmse_deg(i_sep, i_method) = sqrt(mean( ...
            sq_angle_error(:, i_method)));
        channel_nmse_dB(i_sep, i_method) = 10 * log10( ...
            max(mean(nmse_trial(:, i_method)), eps));
        average_gn_iterations(i_sep, i_method) = mean( ...
            iteration_trial(:, i_method));
        average_rejections(i_sep, i_method) = mean( ...
            rejection_trial(:, i_method));
        trace_store{i_sep, i_method} = squeeze( ...
            cost_trace_trial(:, :, i_method));
    end
    initial_doa_rmse_deg(i_sep) = sqrt(mean(initial_angle_sq_error));
    initial_channel_nmse_dB(i_sep) = 10 * log10( ...
        max(mean(initial_nmse_trial), eps));

    fprintf(['GN separation %4.1f deg: init RMSE %.3f deg, ', ...
        'success [%.3f %.3f %.3f], final RMSE [%.3f %.3f %.3f] deg\n'], ...
        separation_deg(i_sep), initial_doa_rmse_deg(i_sep), ...
        success_probability(i_sep, :), doa_rmse_deg(i_sep, :));
end

gn_table = table(separation_deg(:), success_probability(:, 1), ...
    success_probability(:, 2), success_probability(:, 3), ...
    doa_rmse_deg(:, 1), doa_rmse_deg(:, 2), doa_rmse_deg(:, 3), ...
    channel_nmse_dB(:, 1), channel_nmse_dB(:, 2), ...
    channel_nmse_dB(:, 3), average_gn_iterations(:, 1), ...
    average_gn_iterations(:, 2), average_gn_iterations(:, 3), ...
    average_rejections(:, 1), average_rejections(:, 2), ...
    average_rejections(:, 3), initial_doa_rmse_deg, ...
    initial_channel_nmse_dB, ...
    'VariableNames', {'AngularSeparation_deg', ...
    'Undamped_Success', 'LM_Success', 'StabilityAware_Success', ...
    'Undamped_DOA_RMSE_deg', 'LM_DOA_RMSE_deg', ...
    'StabilityAware_DOA_RMSE_deg', 'Undamped_Channel_NMSE_dB', ...
    'LM_Channel_NMSE_dB', 'StabilityAware_Channel_NMSE_dB', ...
    'Undamped_Iterations', 'LM_Iterations', ...
    'StabilityAware_Iterations', 'Undamped_Rejections', ...
    'LM_Rejections', 'StabilityAware_Rejections', ...
    'Algorithm1_Initial_DOA_RMSE_deg', ...
    'Algorithm1_Initial_Channel_NMSE_dB'});
writetable(gn_table, fullfile(result_dir, ...
    'stability_aware_gn_summary.csv'));

fig = publication_figure([80, 80, 920, 350]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
styles = {'^-', 's--', 'o-'};
method_colors = {colors.gray, colors.orange, colors.blue};
nexttile;
for i_method = 1:num_gn_methods
    plot(separation_deg, channel_nmse_dB(:, i_method), ...
        styles{i_method}, 'Color', method_colors{i_method});
    hold on;
end
grid on;
box on;
xlabel('Inter-user angular separation (degree)');
ylabel('Channel NMSE (dB)');
title('(a) Estimation accuracy');
legend(method_names_gn, 'Location', 'southeast');
xticks(separation_deg);

nexttile;
for i_method = 1:num_gn_methods
    plot(separation_deg, success_probability(:, i_method), ...
        styles{i_method}, 'Color', method_colors{i_method});
    hold on;
end
grid on;
box on;
xlabel('Inter-user angular separation (degree)');
ylabel('Successful recovery probability');
title('(b) Basin reliability');
legend(method_names_gn, 'Location', 'southeast');
xticks(separation_deg);
ylim([0, 1.05]);
export_publication_figure(fig, result_dir, ...
    'fig_stability_aware_gn_algorithm1_initialization');

%% 2) Distance-certified finite-alphabet GS under perfect CSI
% This experiment is unaffected by channel initialization and is intentionally
% skipped here.  Its functions remain below because the end-to-end experiment
% uses the same detector.
if false
det_cfg = cfg;
det_cfg.theta_true_deg = [15; 35];
det_cfg.phi_true_deg = [8; -10];
det_cfg.num_trials = 30000;
snr_detection_dB = -15:2:-3;
snr_complexity_dB = -3:3:15;
det_cfg.num_complexity_trials = 10000;

H_data = generate_channel_matrix(det_cfg.theta_true_deg, ...
    det_cfg.phi_true_deg, det_cfg.path_gains, det_cfg);
G_base = H_data;
G_data = repmat(G_base, det_cfg.num_reference_states, 1);
b_data = make_four_phase_data_reference(det_cfg);
mu_codebook = intensity_codebook(G_data, b_data, candidate_symbols);
dmin_exact = minimum_energy_distance(mu_codebook);
dmin_lower = four_phase_distance_lower_bound( ...
    G_base, qpsk, det_cfg.reference_amplitude, det_cfg.num_users);
[ser_union_bound, ser_union_bound_raw] = symbol_error_union_bound( ...
    mu_codebook, candidate_indices, snr_detection_dB);

num_snr_detection = numel(snr_detection_dB);
ser_ml = zeros(num_snr_detection, 1);
ser_linear = zeros(num_snr_detection, 1);
ser_continuous_gs = zeros(num_snr_detection, 1);
ser_certified_gs = zeros(num_snr_detection, 1);
avg_iter_fixed = zeros(num_snr_detection, 1);
avg_iter_exact = zeros(num_snr_detection, 1);
avg_iter_lower = zeros(num_snr_detection, 1);

for i_snr = 1:num_snr_detection
    rng(200000 + i_snr, 'twister');
    num_candidates = size(candidate_symbols, 1);
    true_ids = randi(num_candidates, det_cfg.num_trials, 1);
    noise_variance = mean(mu_codebook(:).^2) / ...
        10^(snr_detection_dB(i_snr) / 10);
    observations = mu_codebook(:, true_ids) + ...
        sqrt(noise_variance) * randn(size(mu_codebook, 1), ...
        det_cfg.num_trials);

    mu_norm = sum(mu_codebook.^2, 1).';
    metrics = mu_norm + sum(observations.^2, 1) - ...
        2 * real(mu_codebook' * observations);
    [~, detected_ml] = min(metrics, [], 1);
    detected_ml = detected_ml(:);

    detected_linear = zeros(det_cfg.num_trials, det_cfg.num_users);
    detected_continuous = zeros(det_cfg.num_trials, det_cfg.num_users);
    detected_certified = zeros(det_cfg.num_trials, det_cfg.num_users);
    [linear_matrix, back_projection] = detector_matrices( ...
        G_data, b_data, det_cfg);

    for trial = 1:det_cfg.num_trials
        z = observations(:, trial);
        s_linear = linearized_symbol_estimate(z, b_data, ...
            linear_matrix, det_cfg.num_users);
        detected_linear(trial, :) = quantize_indices( ...
            s_linear, qpsk).';
        detected_continuous(trial, :) = continuous_affine_gs_detector( ...
            z, G_data, b_data, qpsk, s_linear, ...
            back_projection, det_cfg).';
        [detected_certified(trial, :), ~] = finite_alphabet_gs_detector( ...
            z, G_data, b_data, qpsk, s_linear, ...
            back_projection, dmin_exact, det_cfg, true);
    end

    true_indices = candidate_indices(true_ids, :);
    ml_indices = candidate_indices(detected_ml, :);
    ser_ml(i_snr) = symbol_error_rate(ml_indices, true_indices);
    ser_linear(i_snr) = symbol_error_rate( ...
        detected_linear, true_indices);
    ser_continuous_gs(i_snr) = symbol_error_rate( ...
        detected_continuous, true_indices);
    ser_certified_gs(i_snr) = symbol_error_rate( ...
        detected_certified, true_indices);
    fprintf(['Detection %5.1f dB: ML %.3e, linear %.3e, ', ...
        'continuous GS %.3e, certified GS %.3e\n'], ...
        snr_detection_dB(i_snr), ser_ml(i_snr), ser_linear(i_snr), ...
        ser_continuous_gs(i_snr), ser_certified_gs(i_snr));
end

detection_table = table(snr_detection_dB(:), ser_ml, ser_linear, ...
    ser_continuous_gs, ser_certified_gs, ser_union_bound, ...
    ser_union_bound_raw, ...
    'VariableNames', {'SNR_dB', 'ML_SER', 'Linearized_SER', ...
    'Continuous_GS_SER', 'Certified_FiniteAlphabet_GS_SER', ...
    'Analytical_SER_Union_Bound', 'Unclipped_SER_Union_Bound'});
writetable(detection_table, fullfile(result_dir, ...
    'certified_gs_summary.csv'));

avg_iter_fixed = det_cfg.gs_max_iter * ones(numel(snr_complexity_dB), 1);
avg_iter_exact = zeros(numel(snr_complexity_dB), 1);
avg_iter_lower = zeros(numel(snr_complexity_dB), 1);
for i_snr = 1:numel(snr_complexity_dB)
    rng(210000 + i_snr, 'twister');
    true_ids = randi(size(candidate_symbols, 1), ...
        det_cfg.num_complexity_trials, 1);
    noise_variance = mean(mu_codebook(:).^2) / ...
        10^(snr_complexity_dB(i_snr) / 10);
    observations = mu_codebook(:, true_ids) + ...
        sqrt(noise_variance) * randn(size(mu_codebook, 1), ...
        det_cfg.num_complexity_trials);
    iter_exact = zeros(det_cfg.num_complexity_trials, 1);
    iter_lower = zeros(det_cfg.num_complexity_trials, 1);
    [linear_matrix, back_projection] = detector_matrices( ...
        G_data, b_data, det_cfg);
    for trial = 1:det_cfg.num_complexity_trials
        z = observations(:, trial);
        s_linear = linearized_symbol_estimate(z, b_data, ...
            linear_matrix, det_cfg.num_users);
        [~, iter_exact(trial)] = finite_alphabet_gs_detector( ...
            z, G_data, b_data, qpsk, s_linear, ...
            back_projection, dmin_exact, det_cfg, true);
        [~, iter_lower(trial)] = finite_alphabet_gs_detector( ...
            z, G_data, b_data, qpsk, s_linear, ...
            back_projection, dmin_lower, det_cfg, true);
    end
    avg_iter_exact(i_snr) = mean(iter_exact);
    avg_iter_lower(i_snr) = mean(iter_lower);
    fprintf('Certificate %5.1f dB: iterations %.2f/%.2f/%.2f\n', ...
        snr_complexity_dB(i_snr), avg_iter_fixed(i_snr), ...
        avg_iter_exact(i_snr), avg_iter_lower(i_snr));
end

complexity_table = table(snr_complexity_dB(:), avg_iter_fixed, ...
    avg_iter_exact, avg_iter_lower, ...
    'VariableNames', {'SNR_dB', 'FixedBudget_GS_Iterations', ...
    'ExactCertificate_GS_Iterations', ...
    'FourPhaseBound_GS_Iterations'});
writetable(complexity_table, fullfile(result_dir, ...
    'certified_gs_complexity_summary.csv'));

fig = publication_figure([80, 80, 920, 350]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
semilogy(snr_detection_dB, plotting_ser(ser_ml, ...
    det_cfg.num_trials, det_cfg.num_users), 'o-', 'Color', colors.blue);
hold on;
semilogy(snr_detection_dB, plotting_ser(ser_linear, ...
    det_cfg.num_trials, det_cfg.num_users), '^:', 'Color', colors.gray);
semilogy(snr_detection_dB, plotting_ser(ser_continuous_gs, ...
    det_cfg.num_trials, det_cfg.num_users), 'x-.', 'Color', colors.purple);
semilogy(snr_detection_dB, plotting_ser(ser_certified_gs, ...
    det_cfg.num_trials, det_cfg.num_users), 'd--', 'Color', colors.yellow);
semilogy(snr_detection_dB, max(ser_union_bound, eps), 'k--', ...
    'LineWidth', 1.5);
grid on;
box on;
xlabel('SNR (dB)');
ylabel('SER');
title('(a) Detection accuracy');
legend('ML', 'Linearized detector', 'Continuous affine GS', ...
    'Certified finite-alphabet GS', 'Analytical union bound', ...
    'Location', 'southwest');
ylim([1e-4, 1]);

nexttile;
plot(snr_complexity_dB, avg_iter_fixed, 'x-.', 'Color', colors.purple);
hold on;
plot(snr_complexity_dB, avg_iter_exact, 'd-', 'Color', colors.yellow);
plot(snr_complexity_dB, avg_iter_lower, 's--', 'Color', colors.orange);
grid on;
box on;
xlabel('SNR (dB)');
ylabel('Average GS iterations');
title('(b) Certificate-induced saving');
legend('Fixed iteration budget', 'Exact d_{min}^{(E)} certificate', ...
    'Four-phase lower-bound certificate', 'Location', 'southwest');
ylim([0, det_cfg.gs_max_iter + 1]);
export_publication_figure(fig, result_dir, ...
    'fig_certified_gs_accuracy_complexity');
end

%% 3) End-to-end accuracy and nonlinear workload
e2e_cfg = cfg;
e2e_cfg.num_pilots = 16;
e2e_cfg.theta_true_deg = [18; 22];
e2e_cfg.phi_true_deg = [5; 5];
e2e_cfg.data_snr_dB = 12;
e2e_cfg.num_channel_trials = 100;
e2e_cfg.num_data_trials = 100;
if quick_check
    e2e_cfg.num_channel_trials = 4;
    e2e_cfg.num_data_trials = 10;
end
training_snr_dB = 0:4:20;

[pilot_e2e, reference_training_e2e, pilot_cycles_e2e] = ...
    make_four_phase_training(e2e_cfg, 20260920);
eta_true_e2e = pack_eta(e2e_cfg.theta_true_deg, ...
    e2e_cfg.phi_true_deg, e2e_cfg.path_gains);
[mu_training_e2e, ~, H_true_e2e] = training_model_and_jacobian( ...
    eta_true_e2e, pilot_e2e, reference_training_e2e, e2e_cfg);
G_true_e2e = repmat(H_true_e2e, e2e_cfg.num_reference_states, 1);
b_e2e = make_four_phase_data_reference(e2e_cfg);
mu_true_data = intensity_codebook(G_true_e2e, b_e2e, candidate_symbols);
data_noise_variance = mean(mu_true_data(:).^2) / ...
    10^(e2e_cfg.data_snr_dB / 10);

ser_e2e_baseline = zeros(numel(training_snr_dB), 1);
ser_e2e_proposed = zeros(numel(training_snr_dB), 1);
ser_e2e_perfect_ml = zeros(numel(training_snr_dB), 1);
workload_baseline = zeros(numel(training_snr_dB), 1);
workload_proposed = zeros(numel(training_snr_dB), 1);
channel_nmse_baseline_dB = zeros(numel(training_snr_dB), 1);
channel_nmse_proposed_dB = zeros(numel(training_snr_dB), 1);
initial_nmse_e2e_dB = zeros(numel(training_snr_dB), 1);
initial_doa_rmse_e2e_deg = zeros(numel(training_snr_dB), 1);

for i_snr = 1:numel(training_snr_dB)
    training_noise_variance = mean(mu_training_e2e.^2) / ...
        10^(training_snr_dB(i_snr) / 10);
    error_baseline = 0;
    error_proposed = 0;
    error_perfect = 0;
    iteration_baseline = 0;
    iteration_proposed = 0;
    nmse_baseline = zeros(e2e_cfg.num_channel_trials, 1);
    nmse_proposed = zeros(e2e_cfg.num_channel_trials, 1);
    nmse_initial = zeros(e2e_cfg.num_channel_trials, 1);
    angle_initial = zeros(e2e_cfg.num_channel_trials, 1);

    for channel_trial = 1:e2e_cfg.num_channel_trials
        rng(300000 + 1000 * i_snr + channel_trial, 'twister');
        z_training = mu_training_e2e + sqrt(training_noise_variance) * ...
            randn(size(mu_training_e2e));
        [eta_initial, H_initial] = algorithm1_initialization( ...
            z_training, pilot_cycles_e2e, e2e_cfg, coarse_dictionary, ...
            coarse_theta_deg, coarse_phi_deg);
        nmse_initial(channel_trial) = norm( ...
            H_initial - H_true_e2e, 'fro')^2 / ...
            max(norm(H_true_e2e, 'fro')^2, eps);
        [~, angle_initial(channel_trial)] = matched_doa_error( ...
            eta_initial, e2e_cfg.theta_true_deg, ...
            e2e_cfg.phi_true_deg, e2e_cfg);
        [eta_lm, info_lm] = parametric_gn(z_training, pilot_e2e, ...
            reference_training_e2e, eta_initial, 2, e2e_cfg);
        [eta_stable, info_stable] = parametric_gn(z_training, pilot_e2e, ...
            reference_training_e2e, eta_initial, 3, e2e_cfg);
        H_lm = eta_to_channel(eta_lm, e2e_cfg);
        H_stable = eta_to_channel(eta_stable, e2e_cfg);
        nmse_baseline(channel_trial) = norm(H_lm - H_true_e2e, 'fro')^2 / ...
            max(norm(H_true_e2e, 'fro')^2, eps);
        nmse_proposed(channel_trial) = norm(H_stable - H_true_e2e, 'fro')^2 / ...
            max(norm(H_true_e2e, 'fro')^2, eps);

        G_lm = repmat(H_lm, e2e_cfg.num_reference_states, 1);
        G_stable = repmat(H_stable, e2e_cfg.num_reference_states, 1);
        mu_lm = intensity_codebook(G_lm, b_e2e, candidate_symbols);
        mu_stable = intensity_codebook(G_stable, b_e2e, candidate_symbols);
        dmin_stable = minimum_energy_distance(mu_stable);
        [A_lm, B_lm] = detector_matrices(G_lm, b_e2e, e2e_cfg);
        [A_stable, B_stable] = detector_matrices( ...
            G_stable, b_e2e, e2e_cfg);

        true_ids = randi(size(candidate_symbols, 1), ...
            e2e_cfg.num_data_trials, 1);
        z_data = mu_true_data(:, true_ids) + sqrt(data_noise_variance) * ...
            randn(size(mu_true_data, 1), e2e_cfg.num_data_trials);

        for data_trial = 1:e2e_cfg.num_data_trials
            z = z_data(:, data_trial);
            true_index = candidate_indices(true_ids(data_trial), :);

            linear_lm = linearized_symbol_estimate( ...
                z, b_e2e, A_lm, e2e_cfg.num_users);
            [index_lm, iter_lm] = finite_alphabet_gs_detector( ...
                z, G_lm, b_e2e, qpsk, linear_lm, B_lm, ...
                0, e2e_cfg, false);
            linear_stable = linearized_symbol_estimate( ...
                z, b_e2e, A_stable, e2e_cfg.num_users);
            [index_stable, iter_stable] = finite_alphabet_gs_detector( ...
                z, G_stable, b_e2e, qpsk, linear_stable, B_stable, ...
                dmin_stable, e2e_cfg, true);

            metric_true = sum((mu_true_data - z).^2, 1);
            [~, id_perfect] = min(metric_true);
            index_perfect = candidate_indices(id_perfect, :);
            error_baseline = error_baseline + ...
                sum(index_lm(:).' ~= true_index);
            error_proposed = error_proposed + ...
                sum(index_stable(:).' ~= true_index);
            error_perfect = error_perfect + ...
                sum(index_perfect ~= true_index);
            iteration_baseline = iteration_baseline + iter_lm;
            iteration_proposed = iteration_proposed + iter_stable;
        end
        iteration_baseline = iteration_baseline + info_lm.iterations;
        iteration_proposed = iteration_proposed + info_stable.iterations;
    end

    total_symbol_decisions = e2e_cfg.num_channel_trials * ...
        e2e_cfg.num_data_trials * e2e_cfg.num_users;
    ser_e2e_baseline(i_snr) = error_baseline / total_symbol_decisions;
    ser_e2e_proposed(i_snr) = error_proposed / total_symbol_decisions;
    ser_e2e_perfect_ml(i_snr) = error_perfect / total_symbol_decisions;
    workload_baseline(i_snr) = iteration_baseline / ...
        e2e_cfg.num_channel_trials;
    workload_proposed(i_snr) = iteration_proposed / ...
        e2e_cfg.num_channel_trials;
    channel_nmse_baseline_dB(i_snr) = 10 * log10( ...
        max(mean(nmse_baseline), eps));
    channel_nmse_proposed_dB(i_snr) = 10 * log10( ...
        max(mean(nmse_proposed), eps));
    initial_nmse_e2e_dB(i_snr) = 10 * log10( ...
        max(mean(nmse_initial), eps));
    initial_doa_rmse_e2e_deg(i_snr) = sqrt(mean(angle_initial));

    fprintf(['End-to-end training %4.1f dB: SER baseline %.3e, ', ...
        'proposed %.3e, perfect ML %.3e; init RMSE %.2f deg\n'], ...
        training_snr_dB(i_snr), ser_e2e_baseline(i_snr), ...
        ser_e2e_proposed(i_snr), ser_e2e_perfect_ml(i_snr), ...
        initial_doa_rmse_e2e_deg(i_snr));
end

e2e_table = table(training_snr_dB(:), ser_e2e_baseline, ...
    ser_e2e_proposed, ser_e2e_perfect_ml, workload_baseline, ...
    workload_proposed, channel_nmse_baseline_dB, ...
    channel_nmse_proposed_dB, initial_doa_rmse_e2e_deg, ...
    initial_nmse_e2e_dB, ...
    'VariableNames', {'Training_SNR_dB', 'Conventional_EndToEnd_SER', ...
    'TheoryGuided_EndToEnd_SER', 'PerfectCSI_ML_SER', ...
    'Conventional_Nonlinear_Iterations_perBlock', ...
    'TheoryGuided_Nonlinear_Iterations_perBlock', ...
    'Conventional_Channel_NMSE_dB', ...
    'TheoryGuided_Channel_NMSE_dB', ...
    'Algorithm1_Initial_DOA_RMSE_deg', ...
    'Algorithm1_Initial_Channel_NMSE_dB'});
writetable(e2e_table, fullfile(result_dir, ...
    'theory_guided_end_to_end_summary.csv'));

fig = publication_figure([80, 80, 920, 350]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
semilogy(training_snr_dB, plotting_ser(ser_e2e_baseline, ...
    e2e_cfg.num_channel_trials * e2e_cfg.num_data_trials, ...
    e2e_cfg.num_users), 's--', 'Color', colors.orange);
hold on;
semilogy(training_snr_dB, plotting_ser(ser_e2e_proposed, ...
    e2e_cfg.num_channel_trials * e2e_cfg.num_data_trials, ...
    e2e_cfg.num_users), 'o-', 'Color', colors.blue);
semilogy(training_snr_dB, plotting_ser(ser_e2e_perfect_ml, ...
    e2e_cfg.num_channel_trials * e2e_cfg.num_data_trials, ...
    e2e_cfg.num_users), 'k-.');
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
plot(training_snr_dB, channel_nmse_baseline_dB, ...
    's--', 'Color', colors.orange);
hold on;
plot(training_snr_dB, channel_nmse_proposed_dB, ...
    'o-', 'Color', colors.blue);
grid on;
box on;
xlabel('Training SNR (dB)');
ylabel('Channel NMSE (dB)');
title('(b) Channel-estimation accuracy');
legend('Conventional LM--GN', 'Stability-aware GN', ...
    'Location', 'southwest');
export_publication_figure(fig, result_dir, ...
    'fig_theory_guided_end_to_end_algorithm1_initialization');

%% Diagnostic comparison with the current oracle-initialized results
old_result_dir = fullfile(fileparts(script_dir), ...
    'ch5_theory_guided_algorithms_2026', 'results');
old_gn = readtable(fullfile(old_result_dir, ...
    'stability_aware_gn_summary.csv'));
old_e2e = readtable(fullfile(old_result_dir, ...
    'theory_guided_end_to_end_summary.csv'));

gn_comparison = table(separation_deg(:), ...
    old_gn.StabilityAware_Success, success_probability(:, 3), ...
    old_gn.StabilityAware_Channel_NMSE_dB, channel_nmse_dB(:, 3), ...
    initial_doa_rmse_deg, initial_channel_nmse_dB, ...
    'VariableNames', {'AngularSeparation_deg', ...
    'OldOracle_RankOnePilot_Success', ...
    'Algorithm1_FullRankPilot_Success', ...
    'OldOracle_RankOnePilot_Channel_NMSE_dB', ...
    'Algorithm1_FullRankPilot_Channel_NMSE_dB', ...
    'Algorithm1_Initial_DOA_RMSE_deg', ...
    'Algorithm1_Initial_Channel_NMSE_dB'});
writetable(gn_comparison, fullfile(result_dir, ...
    'comparison_old_oracle_vs_algorithm1_gn.csv'));

e2e_comparison = table(training_snr_dB(:), ...
    old_e2e.TheoryGuided_EndToEnd_SER, ser_e2e_proposed, ...
    old_e2e.TheoryGuided_Channel_NMSE_dB, ...
    channel_nmse_proposed_dB, initial_doa_rmse_e2e_deg, ...
    'VariableNames', {'Training_SNR_dB', ...
    'OldOracle_RankOnePilot_SER', 'Algorithm1_FullRankPilot_SER', ...
    'OldOracle_RankOnePilot_Channel_NMSE_dB', ...
    'Algorithm1_FullRankPilot_Channel_NMSE_dB', ...
    'Algorithm1_Initial_DOA_RMSE_deg'});
writetable(e2e_comparison, fullfile(result_dir, ...
    'comparison_old_oracle_vs_algorithm1_end_to_end.csv'));

fig = publication_figure([80, 80, 920, 350]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(separation_deg, old_gn.StabilityAware_Channel_NMSE_dB, ...
    's--', 'Color', colors.orange);
hold on;
plot(separation_deg, channel_nmse_dB(:, 3), ...
    'o-', 'Color', colors.blue);
grid on;
box on;
xlabel('Inter-user angular separation (degree)');
ylabel('Channel NMSE (dB)');
title('(a) Stability-aware GN accuracy');
legend('Old oracle / rank-one pilots', ...
    'Algorithm 1 / full-rank pilots', 'Location', 'southeast');
xticks(separation_deg);
nexttile;
plot(separation_deg, old_gn.StabilityAware_Success, ...
    's--', 'Color', colors.orange);
hold on;
plot(separation_deg, success_probability(:, 3), ...
    'o-', 'Color', colors.blue);
grid on;
box on;
xlabel('Inter-user angular separation (degree)');
ylabel('Successful recovery probability');
title('(b) Basin reliability');
legend('Old oracle / rank-one pilots', ...
    'Algorithm 1 / full-rank pilots', 'Location', 'southeast');
xticks(separation_deg);
ylim([0, 1.05]);
export_publication_figure(fig, result_dir, ...
    'fig_diagnostic_old_oracle_vs_algorithm1_gn');

fig = publication_figure([80, 80, 920, 350]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
semilogy(training_snr_dB, plotting_ser( ...
    old_e2e.TheoryGuided_EndToEnd_SER, ...
    e2e_cfg.num_channel_trials * e2e_cfg.num_data_trials, ...
    e2e_cfg.num_users), 's--', 'Color', colors.orange);
hold on;
semilogy(training_snr_dB, plotting_ser(ser_e2e_proposed, ...
    e2e_cfg.num_channel_trials * e2e_cfg.num_data_trials, ...
    e2e_cfg.num_users), 'o-', 'Color', colors.blue);
grid on;
box on;
xlabel('Training SNR (dB)');
ylabel('End-to-end SER');
title('(a) Theory-guided receiver');
legend('Old oracle / rank-one pilots', ...
    'Algorithm 1 / full-rank pilots', 'Location', 'southwest');
nexttile;
plot(training_snr_dB, old_e2e.TheoryGuided_Channel_NMSE_dB, ...
    's--', 'Color', colors.orange);
hold on;
plot(training_snr_dB, channel_nmse_proposed_dB, ...
    'o-', 'Color', colors.blue);
grid on;
box on;
xlabel('Training SNR (dB)');
ylabel('Channel NMSE (dB)');
title('(b) Channel-estimation accuracy');
legend('Old oracle / rank-one pilots', ...
    'Algorithm 1 / full-rank pilots', 'Location', 'southwest');
export_publication_figure(fig, result_dir, ...
    'fig_diagnostic_old_oracle_vs_algorithm1_end_to_end');

save(fullfile(result_dir, 'algorithm1_practical_initialization_results.mat'), ...
    'cfg', 'gn_cfg', 'e2e_cfg', 'separation_deg', ...
    'success_probability', 'doa_rmse_deg', 'channel_nmse_dB', ...
    'average_gn_iterations', 'average_rejections', ...
    'initial_doa_rmse_deg', 'initial_channel_nmse_dB', ...
    'pilot_rank', 'pilot_condition', 'training_snr_dB', ...
    'ser_e2e_baseline', 'ser_e2e_proposed', 'ser_e2e_perfect_ml', ...
    'workload_baseline', 'workload_proposed', ...
    'channel_nmse_baseline_dB', 'channel_nmse_proposed_dB', ...
    'initial_doa_rmse_e2e_deg', 'initial_nmse_e2e_dB');

fprintf('\nAll outputs saved to:\n%s\n', result_dir);

%% Local functions
function [S, B, pilot_cycles] = make_four_phase_training(cfg, seed)
    num_states = cfg.num_reference_states;
    assert(num_states == 4, ...
        'The explicit field initializer requires four phase states.');
    assert(mod(cfg.num_pilots, num_states) == 0, ...
        'T_p must be divisible by four.');
    num_cycles = cfg.num_pilots / num_states;
    pilot_cycles = make_correlated_full_rank_pilots( ...
        num_cycles, cfg.num_users, cfg.pilot_shared_fraction, seed);
    S = repelem(pilot_cycles, num_states, 1);
    reference_phase = repmat((0:num_states - 1).' * pi / 2, ...
        num_cycles, 1);
    space_phase = 2 * pi * (0:cfg.num_ris - 1) / cfg.num_ris;
    B = cfg.reference_amplitude * exp( ...
        1j * (reference_phase + space_phase));
end

function S = make_correlated_full_rank_pilots(num_cycles, users, ...
        shared_fraction, seed)
    old_rng = rng;
    rng(seed, 'twister');
    qpsk = [1 + 1j; -1 + 1j; -1 - 1j; 1 - 1j] / sqrt(2);
    base = qpsk(randi(4, num_cycles, 1));
    S = repmat(base, 1, users);
    for user = 2:users
        num_changes = max(1, round((1 - shared_fraction) * num_cycles));
        changed = randperm(num_cycles, num_changes);
        S(changed, user) = 1j * S(changed, user);
    end
    assert(rank(S) == users, ...
        'The user-pilot matrix must have full column rank.');
    rng(old_rng);
end

function [dictionary, theta_deg, phi_deg] = make_coarse_dictionary(cfg)
    [theta_mesh, phi_mesh] = ndgrid( ...
        cfg.coarse_theta_grid_deg, cfg.coarse_phi_grid_deg);
    theta_deg = theta_mesh(:);
    phi_deg = phi_mesh(:);
    dictionary = zeros(cfg.num_ris, numel(theta_deg));
    for grid_index = 1:numel(theta_deg)
        dictionary(:, grid_index) = steering_vector_with_derivatives( ...
            theta_deg(grid_index), phi_deg(grid_index), cfg);
    end
end

function [eta_initial, H_initial] = algorithm1_initialization( ...
        z, pilot_cycles, cfg, dictionary, theta_grid_deg, phi_grid_deg)
    num_cycles = size(pilot_cycles, 1);
    Z = reshape(z, cfg.num_pilots, cfg.num_ris);
    space_phase = 2 * pi * (0:cfg.num_ris - 1) / cfg.num_ris;
    field_cycles = zeros(num_cycles, cfg.num_ris);
    for cycle = 1:num_cycles
        rows = (cycle - 1) * 4 + (1:4);
        rotated_field = (Z(rows(1), :) - Z(rows(3), :) + ...
            1j * (Z(rows(2), :) - Z(rows(4), :))) / ...
            (4 * cfg.reference_amplitude);
        field_cycles(cycle, :) = rotated_field .* exp(1j * space_phase);
    end

    H_unstructured = (pilot_cycles \ field_cycles).';
    theta_initial = zeros(cfg.num_users, 1);
    phi_initial = zeros(cfg.num_users, 1);
    alpha_initial = zeros(cfg.num_users, 1);
    dictionary_norm = sum(abs(dictionary).^2, 1).';
    for user = 1:cfg.num_users
        score = abs(dictionary' * H_unstructured(:, user)).^2 ./ ...
            max(dictionary_norm, eps);
        [~, best] = max(score);
        theta_initial(user) = theta_grid_deg(best);
        phi_initial(user) = phi_grid_deg(best);
        steering = dictionary(:, best);
        alpha_initial(user) = (steering' * H_unstructured(:, user)) / ...
            (steering' * steering);
    end
    eta_initial = project_eta(pack_eta( ...
        theta_initial, phi_initial, alpha_initial), cfg);
    H_initial = eta_to_channel(eta_initial, cfg);
end

function b = make_four_phase_data_reference(cfg)
    phases = (0:cfg.num_reference_states - 1).' * ...
        2 * pi / cfg.num_reference_states;
    B = cfg.reference_amplitude * exp(1j * phases) * ...
        ones(1, cfg.num_ris);
    b = reshape(B.', [], 1);
end

function eta = pack_eta(theta_deg, phi_deg, alpha)
    eta = [theta_deg(:) * pi / 180; phi_deg(:) * pi / 180; ...
        real(alpha(:)); imag(alpha(:))];
end

function eta = project_eta(eta, cfg)
    L = cfg.num_users;
    eta(1:L) = min(max(eta(1:L), 0), 60 * pi / 180);
    eta(L + 1:2 * L) = min(max(eta(L + 1:2 * L), ...
        -30 * pi / 180), 30 * pi / 180);
    gain = eta(2 * L + 1:3 * L) + 1j * eta(3 * L + 1:4 * L);
    magnitude = abs(gain);
    too_large = magnitude > 2;
    gain(too_large) = 2 * gain(too_large) ./ magnitude(too_large);
    eta(2 * L + 1:3 * L) = real(gain);
    eta(3 * L + 1:4 * L) = imag(gain);
end

function H = eta_to_channel(eta, cfg)
    L = cfg.num_users;
    theta_deg = eta(1:L) * 180 / pi;
    phi_deg = eta(L + 1:2 * L) * 180 / pi;
    alpha = eta(2 * L + 1:3 * L) + ...
        1j * eta(3 * L + 1:4 * L);
    H = generate_channel_matrix(theta_deg, phi_deg, alpha, cfg);
end

function H = generate_channel_matrix(theta_deg, phi_deg, alpha, cfg)
    H = zeros(cfg.num_ris, cfg.num_users);
    for user = 1:cfg.num_users
        v = steering_vector_with_derivatives( ...
            theta_deg(user), phi_deg(user), cfg);
        H(:, user) = alpha(user) * v;
    end
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

function [mu, J, H] = training_model_and_jacobian(eta, S, B, cfg)
    L = cfg.num_users;
    theta_deg = eta(1:L) * 180 / pi;
    phi_deg = eta(L + 1:2 * L) * 180 / pi;
    alpha = eta(2 * L + 1:3 * L) + ...
        1j * eta(3 * L + 1:4 * L);
    H = zeros(cfg.num_ris, L);
    field_derivative = zeros(numel(B), 4 * L);
    field = zeros(size(B));

    for user = 1:L
        [v, dv_theta, dv_phi] = steering_vector_with_derivatives( ...
            theta_deg(user), phi_deg(user), cfg);
        H(:, user) = alpha(user) * v;
        field = field + S(:, user) * H(:, user).';
        derivatives = {alpha(user) * dv_theta, ...
            alpha(user) * dv_phi, v, 1j * v};
        columns = [user, L + user, 2 * L + user, 3 * L + user];
        for derivative_index = 1:4
            derivative_matrix = S(:, user) * derivatives{derivative_index}.';
            field_derivative(:, columns(derivative_index)) = ...
                derivative_matrix(:);
        end
    end

    affine_field = field + B;
    mu = abs(affine_field(:)).^2;
    J = 2 * real(conj(affine_field(:)) .* field_derivative);
end

function [eta, information] = parametric_gn(z, S, B, eta_initial, ...
        method, cfg)
    eta = eta_initial;
    lambda_lm = cfg.gn_initial_damping;
    [mu, J] = training_model_and_jacobian(eta, S, B, cfg);
    residual = z - mu;
    cost = 0.5 * norm(residual)^2;
    cost_trace = nan(1, cfg.gn_max_iter + 1);
    cost_trace(1) = cost;
    rejections = 0;

    for iteration = 1:cfg.gn_max_iter
        singular_values = svd(J, 'econ');
        sigma_max = max(singular_values);
        sigma_min = min(singular_values);
        if method == 1
            lambda = 0;
        elseif method == 2
            lambda = lambda_lm;
        else
            lambda_stability = max((sigma_max^2 - ...
                cfg.gn_target_condition * sigma_min^2) / ...
                (cfg.gn_target_condition - 1), 0);
            lambda = max(lambda_lm, lambda_stability);
        end

        increment = (J.' * J + (lambda + 1e-10) * eye(size(J, 2))) \ ...
            (J.' * residual);
        eta_trial = project_eta(eta + increment, cfg);
        [mu_trial, J_trial] = training_model_and_jacobian( ...
            eta_trial, S, B, cfg);
        residual_trial = z - mu_trial;
        cost_trial = 0.5 * norm(residual_trial)^2;

        if method == 1
            eta = eta_trial;
            mu = mu_trial;
            J = J_trial;
            residual = residual_trial;
            cost = cost_trial;
        elseif cost_trial < cost
            eta = eta_trial;
            mu = mu_trial;
            J = J_trial;
            residual = residual_trial;
            cost = cost_trial;
            lambda_lm = max(lambda_lm / 3, 1e-10);
        else
            lambda_lm = min(lambda_lm * 10, 1e10);
            rejections = rejections + 1;
        end

        cost_trace(iteration + 1) = cost;
        if norm(increment) <= cfg.gn_tolerance * (norm(eta) + 1)
            break;
        end
    end

    cost_trace(iteration + 2:end) = cost;
    information.iterations = iteration;
    information.rejections = rejections;
    information.cost_trace = cost_trace;
end

function [rmse_deg, squared_error] = matched_doa_error( ...
        eta, theta_true_deg, phi_true_deg, cfg)
    L = cfg.num_users;
    theta_hat = eta(1:L) * 180 / pi;
    phi_hat = eta(L + 1:2 * L) * 180 / pi;
    direct = mean((theta_hat - theta_true_deg).^2 + ...
        (phi_hat - phi_true_deg).^2);
    swapped = mean((flipud(theta_hat) - theta_true_deg).^2 + ...
        (flipud(phi_hat) - phi_true_deg).^2);
    squared_error = min(direct, swapped);
    rmse_deg = sqrt(squared_error);
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

function lower_bound = four_phase_distance_lower_bound( ...
        G, constellation, reference_amplitude, users)
    [~, candidates] = enumerate_symbol_vectors(constellation, users);
    minimum_field_distance = inf;
    for first = 1:size(candidates, 1)
        for second = first + 1:size(candidates, 1)
            field_distance = norm(G * ...
                (candidates(first, :) - candidates(second, :)).');
            minimum_field_distance = min( ...
                minimum_field_distance, field_distance);
        end
    end
    lower_bound = 2 * sqrt(2) * reference_amplitude * ...
        minimum_field_distance;
end

function [A, back_projection] = detector_matrices(G, b, cfg)
    users = size(G, 2);
    weighted = conj(b) .* G;
    A = 2 * [real(weighted), -imag(weighted)];
    back_projection = (G' * G + cfg.gs_regularization * eye(users)) \ G';
end

function estimate = linearized_symbol_estimate(z, b, A, users)
    regularization = 1e-3 * trace(A.' * A) / max(2 * users, 1);
    real_estimate = (A.' * A + regularization * eye(2 * users)) \ ...
        (A.' * (z - abs(b).^2));
    estimate = real_estimate(1:users) + ...
        1j * real_estimate(users + 1:end);
end

function detected_indices = continuous_affine_gs_detector( ...
        z, G, b, constellation, initial, back_projection, cfg)
    s = initial;
    measured_magnitude = sqrt(max(z, 0));
    for iteration = 1:cfg.gs_max_iter
        affine_field = G * s + b;
        phase = affine_field ./ max(abs(affine_field), eps);
        s_next = back_projection * (measured_magnitude .* phase - b);
        if norm(s_next - s) / max(norm(s), 1) < cfg.gs_tolerance
            s = s_next;
            break;
        end
        s = s_next;
    end
    detected_indices = quantize_indices(s, constellation);
end

function [detected_indices, iterations] = finite_alphabet_gs_detector( ...
        z, G, b, constellation, initial, back_projection, ...
        certificate_distance, cfg, use_certificate)
    indices = quantize_indices(initial, constellation);
    s = constellation(indices);
    measured_magnitude = sqrt(max(z, 0));
    iterations = 0;

    for iteration = 1:cfg.gs_max_iter
        if use_certificate
            residual = norm(z - abs(G * s + b).^2);
            if residual < certificate_distance / 2
                iterations = iteration - 1;
                detected_indices = indices;
                return;
            end
        end
        affine_field = G * s + b;
        phase = affine_field ./ max(abs(affine_field), eps);
        continuous_update = back_projection * ...
            (measured_magnitude .* phase - b);
        indices_next = quantize_indices(continuous_update, constellation);
        s_next = constellation(indices_next);
        s = s_next;
        indices = indices_next;
        iterations = iteration;
    end
    detected_indices = indices;
end

function indices = quantize_indices(symbols, constellation)
    symbols = symbols(:);
    constellation = constellation(:).';
    [~, indices] = min(abs(symbols - constellation), [], 2);
end

function ser = symbol_error_rate(detected, truth)
    ser = mean(detected(:) ~= truth(:));
end

function [bound, raw_bound] = symbol_error_union_bound( ...
        mu, combination_indices, snr_dB)
    num_candidates = size(mu, 2);
    num_symbols = size(combination_indices, 2);
    squared_norms = sum(mu.^2, 1);
    distance = sqrt(max(squared_norms.' + squared_norms - ...
        2 * real(mu' * mu), 0));
    hamming = zeros(num_candidates, num_candidates);
    for candidate = 1:num_candidates
        hamming(candidate, :) = sum(combination_indices ~= ...
            combination_indices(candidate, :), 2).';
    end
    rms_intensity_squared = mean(mu(:).^2);
    raw_bound = zeros(numel(snr_dB), 1);
    for i_snr = 1:numel(snr_dB)
        noise_variance = rms_intensity_squared / 10^(snr_dB(i_snr) / 10);
        pairwise_probability = 0.5 * erfc( ...
            distance / (2 * sqrt(2 * noise_variance)));
        pairwise_probability(1:num_candidates + 1:end) = 0;
        raw_bound(i_snr) = sum(hamming .* pairwise_probability, 'all') / ...
            (num_candidates * num_symbols);
    end
    bound = min(raw_bound, 1);
end

function values = plotting_ser(values, num_trials, num_users)
    values = max(values, 0.5 / (num_trials * num_users));
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
