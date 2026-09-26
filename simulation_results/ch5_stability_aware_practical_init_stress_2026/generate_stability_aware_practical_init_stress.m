%% Stability-aware GN stress test with Algorithm-1 initialization
clc;
clear;
close all;

script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

rng(20260924, 'twister');
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');
quick_check = strcmp(getenv('HOLO_QUICK_CHECK'), '1');

%% Configuration
cfg.num_users = 2;
cfg.ris_rows = 8;
cfg.ris_cols = 8;
cfg.num_ris = cfg.ris_rows * cfg.ris_cols;
cfg.lambda = 1;
cfg.spacing = cfg.lambda / 2;
cfg.reference_amplitude = 1.5;
cfg.num_reference_states = 4;
cfg.num_pilots = 32;
cfg.num_cycles = cfg.num_pilots / cfg.num_reference_states;
cfg.path_gains = [exp(1j * 0.35); 0.85 * exp(-1j * 0.55)];
cfg.theta_true_deg = [19; 21];
cfg.phi_true_deg = [5; 5];
cfg.snr_dB = 8;
cfg.success_threshold_deg = 1;
cfg.gn_max_iter = 40;
cfg.gn_tolerance = 1e-6;
cfg.gn_initial_damping = 1e-3;
cfg.gn_target_condition = 1000;
cfg.coarse_theta_grid_deg = 0:0.5:60;
cfg.coarse_phi_grid_deg = -30:1:30;
cfg.num_monte_carlo = 300;
if quick_check
    cfg.num_monte_carlo = 30;
end

pilot_correlation = [0, 0.5, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995];
method_names = {'Conventional LM', 'Stability-aware GN'};
num_methods = numel(method_names);
num_cases = numel(pilot_correlation);

colors.orange = [0.8500, 0.3250, 0.0980];
colors.blue = [0.0000, 0.4470, 0.7410];
colors.gray = [0.35, 0.35, 0.35];

[coarse_dictionary, coarse_theta_deg, coarse_phi_deg] = ...
    make_coarse_dictionary(cfg);
eta_true = pack_eta(cfg.theta_true_deg, cfg.phi_true_deg, cfg.path_gains);
H_true = eta_to_channel(eta_true, cfg);

pilot_gram_condition = zeros(num_cases, 1);
median_initial_normal_condition = zeros(num_cases, 1);
initial_doa_rmse_deg = zeros(num_cases, 1);
initial_channel_nmse_dB = zeros(num_cases, 1);
success_probability = zeros(num_cases, num_methods);
doa_rmse_deg = zeros(num_cases, num_methods);
channel_nmse_dB = zeros(num_cases, num_methods);
catastrophic_failure_probability = zeros(num_cases, num_methods);
average_iterations = zeros(num_cases, num_methods);
average_rejections = zeros(num_cases, num_methods);
stability_activation_probability = zeros(num_cases, 1);

fprintf('\nStability-aware GN with practical initialization\n');
fprintf('T_p=%d, SNR=%.1f dB, DOA separation=%.1f degree, trials=%d\n', ...
    cfg.num_pilots, cfg.snr_dB, diff(cfg.theta_true_deg), ...
    cfg.num_monte_carlo);

for case_index = 1:num_cases
    rho_p = pilot_correlation(case_index);
    [pilot_symbols, training_reference, pilot_cycles] = ...
        make_four_phase_training(cfg, rho_p, 20261001);
    pilot_gram_condition(case_index) = cond(pilot_cycles' * pilot_cycles);
    [mu_true, ~] = training_model_and_jacobian( ...
        eta_true, pilot_symbols, training_reference, cfg);
    noise_variance = mean(mu_true.^2) / 10^(cfg.snr_dB / 10);

    initial_angle_sq_error = zeros(cfg.num_monte_carlo, 1);
    initial_nmse = zeros(cfg.num_monte_carlo, 1);
    initial_normal_condition = zeros(cfg.num_monte_carlo, 1);
    final_angle_sq_error = zeros(cfg.num_monte_carlo, num_methods);
    final_nmse = zeros(cfg.num_monte_carlo, num_methods);
    iteration_count = zeros(cfg.num_monte_carlo, num_methods);
    rejection_count = zeros(cfg.num_monte_carlo, num_methods);
    activation_count = zeros(cfg.num_monte_carlo, 1);

    for trial = 1:cfg.num_monte_carlo
        rng(1000000 + trial, 'twister');
        z = mu_true + sqrt(noise_variance) * randn(size(mu_true));
        [eta_initial, H_initial] = algorithm1_initialization( ...
            z, pilot_cycles, cfg, coarse_dictionary, ...
            coarse_theta_deg, coarse_phi_deg);
        [~, initial_angle_sq_error(trial)] = matched_doa_error( ...
            eta_initial, cfg.theta_true_deg, cfg.phi_true_deg, cfg);
        initial_nmse(trial) = norm(H_initial - H_true, 'fro')^2 / ...
            max(norm(H_true, 'fro')^2, eps);
        [~, J_initial] = training_model_and_jacobian( ...
            eta_initial, pilot_symbols, training_reference, cfg);
        singular_values = svd(J_initial, 'econ');
        initial_normal_condition(trial) = (max(singular_values) / ...
            max(min(singular_values), eps))^2;

        for method = 1:num_methods
            [eta_hat, information] = parametric_gn( ...
                z, pilot_symbols, training_reference, eta_initial, ...
                method, cfg);
            H_hat = eta_to_channel(eta_hat, cfg);
            [~, final_angle_sq_error(trial, method)] = matched_doa_error( ...
                eta_hat, cfg.theta_true_deg, cfg.phi_true_deg, cfg);
            final_nmse(trial, method) = norm(H_hat - H_true, 'fro')^2 / ...
                max(norm(H_true, 'fro')^2, eps);
            iteration_count(trial, method) = information.iterations;
            rejection_count(trial, method) = information.rejections;
            if method == 2
                activation_count(trial) = information.stability_activations > 0;
            end
        end
    end

    initial_doa_rmse_deg(case_index) = sqrt(mean(initial_angle_sq_error));
    initial_channel_nmse_dB(case_index) = 10 * log10( ...
        max(mean(initial_nmse), eps));
    median_initial_normal_condition(case_index) = median( ...
        initial_normal_condition);
    stability_activation_probability(case_index) = mean(activation_count);
    for method = 1:num_methods
        per_trial_rmse = sqrt(final_angle_sq_error(:, method));
        success_probability(case_index, method) = mean( ...
            per_trial_rmse < cfg.success_threshold_deg);
        catastrophic_failure_probability(case_index, method) = mean( ...
            per_trial_rmse > 5);
        doa_rmse_deg(case_index, method) = sqrt(mean( ...
            final_angle_sq_error(:, method)));
        channel_nmse_dB(case_index, method) = 10 * log10( ...
            max(mean(final_nmse(:, method)), eps));
        average_iterations(case_index, method) = mean( ...
            iteration_count(:, method));
        average_rejections(case_index, method) = mean( ...
            rejection_count(:, method));
    end

    fprintf(['rho %.3f, cond %.1f, median kappa(J''J) %.2e: ', ...
        'success LM/stable %.3f/%.3f, NMSE %.2f/%.2f dB, ', ...
        'activation %.3f\n'], rho_p, pilot_gram_condition(case_index), ...
        median_initial_normal_condition(case_index), ...
        success_probability(case_index, :), ...
        channel_nmse_dB(case_index, :), ...
        stability_activation_probability(case_index));
end

%% Save numerical results
summary_table = table(pilot_correlation(:), pilot_gram_condition, ...
    median_initial_normal_condition, initial_doa_rmse_deg, ...
    initial_channel_nmse_dB, success_probability(:, 1), ...
    success_probability(:, 2), catastrophic_failure_probability(:, 1), ...
    catastrophic_failure_probability(:, 2), doa_rmse_deg(:, 1), ...
    doa_rmse_deg(:, 2), channel_nmse_dB(:, 1), ...
    channel_nmse_dB(:, 2), average_iterations(:, 1), ...
    average_iterations(:, 2), average_rejections(:, 1), ...
    average_rejections(:, 2), stability_activation_probability, ...
    'VariableNames', {'PilotCorrelation', 'PilotGramCondition', ...
    'MedianInitialNormalCondition', 'InitialDOA_RMSE_deg', ...
    'InitialChannel_NMSE_dB', 'LM_Success', ...
    'StabilityAware_Success', 'LM_CatastrophicFailure', ...
    'StabilityAware_CatastrophicFailure', 'LM_DOA_RMSE_deg', ...
    'StabilityAware_DOA_RMSE_deg', 'LM_Channel_NMSE_dB', ...
    'StabilityAware_Channel_NMSE_dB', 'LM_Iterations', ...
    'StabilityAware_Iterations', 'LM_Rejections', ...
    'StabilityAware_Rejections', 'StabilityActivationProbability'});
writetable(summary_table, fullfile(result_dir, ...
    'stability_aware_practical_init_stress_summary.csv'));

%% Publication figures
fig = publication_figure([80, 80, 920, 350]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
semilogx(pilot_gram_condition, success_probability(:, 1), ...
    's--', 'Color', colors.orange);
hold on;
semilogx(pilot_gram_condition, success_probability(:, 2), ...
    'o-', 'Color', colors.blue);
grid on;
box on;
xlabel('Pilot Gram-matrix condition number');
ylabel('Successful recovery probability');
title('(a) Basin reliability');
legend(method_names, 'Location', 'southwest');
ylim([0.65, 1.01]);
nexttile;
semilogx(pilot_gram_condition, channel_nmse_dB(:, 1), ...
    's--', 'Color', colors.orange);
hold on;
semilogx(pilot_gram_condition, channel_nmse_dB(:, 2), ...
    'o-', 'Color', colors.blue);
grid on;
box on;
xlabel('Pilot Gram-matrix condition number');
ylabel('Channel NMSE (dB)');
title('(b) Estimation accuracy');
legend(method_names, 'Location', 'northwest');
export_publication_figure(fig, result_dir, ...
    'fig_stability_benefit_with_algorithm1_initialization');

fig = publication_figure([80, 80, 920, 350]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
loglog(pilot_gram_condition, median_initial_normal_condition, ...
    'd-', 'Color', colors.gray);
hold on;
yline(cfg.gn_target_condition, 'k--', '\kappa_{tar}');
grid on;
box on;
xlabel('Pilot Gram-matrix condition number');
ylabel('Median \kappa(J_0^T J_0)');
title('(a) Initial local conditioning');
nexttile;
semilogx(pilot_gram_condition, average_rejections(:, 1), ...
    's--', 'Color', colors.orange);
hold on;
semilogx(pilot_gram_condition, average_rejections(:, 2), ...
    'o-', 'Color', colors.blue);
grid on;
box on;
xlabel('Pilot Gram-matrix condition number');
ylabel('Average rejected trial steps');
title('(b) Suppression of unreliable updates');
legend(method_names, 'Location', 'northwest');
export_publication_figure(fig, result_dir, ...
    'fig_initial_condition_and_stability_activation');

%% Stability-guided adaptive training
% Start from four pilot cycles (16 intensity snapshots).  When the Jacobian
% at the Algorithm-1 initial point is too ill-conditioned, acquire four
% complementary pilot cycles whose cross-correlation cancels that of the
% initial block.  Conventional LM is used after every initialization so that
% the gain below is attributable to the stability-guided training decision.
adaptive_cfg = cfg;
adaptive_cfg.num_cycles = 4;
adaptive_cfg.num_pilots = 16;
full_cfg = cfg;
full_cfg.num_cycles = 8;
full_cfg.num_pilots = 32;
adaptive_jacobian_trigger = 1000;
adaptive_pilot_trigger = 100;
adaptive_correlation = pilot_correlation;
num_adaptive_cases = numel(adaptive_correlation);
adaptive_pilot_condition = zeros(num_adaptive_cases, 1);
adaptive_median_initial_condition = zeros(num_adaptive_cases, 1);
extension_probability = zeros(num_adaptive_cases, 1);
average_training_snapshots = zeros(num_adaptive_cases, 1);
adaptive_success = zeros(num_adaptive_cases, 3);
adaptive_nmse_dB = zeros(num_adaptive_cases, 3);
adaptive_method_names = {'Fixed T_p=16', ...
    'Stability-guided adaptive', 'Fixed T_p=32'};

for case_index = 1:num_adaptive_cases
    rho_p = adaptive_correlation(case_index);
    [initial_cycles, complementary_cycles] = ...
        make_complementary_pilot_cycles( ...
        adaptive_cfg, rho_p, 20261011);
    full_cycles = [initial_cycles; complementary_cycles];
    [S_initial, B_initial] = build_four_phase_training( ...
        initial_cycles, adaptive_cfg);
    [S_extra, B_extra] = build_four_phase_training( ...
        complementary_cycles, adaptive_cfg);
    [S_full, B_full] = build_four_phase_training(full_cycles, full_cfg);
    adaptive_pilot_condition(case_index) = cond( ...
        initial_cycles' * initial_cycles);
    [mu_initial, ~] = training_model_and_jacobian( ...
        eta_true, S_initial, B_initial, adaptive_cfg);
    [mu_extra, ~] = training_model_and_jacobian( ...
        eta_true, S_extra, B_extra, adaptive_cfg);
    [mu_full, ~] = training_model_and_jacobian( ...
        eta_true, S_full, B_full, full_cfg);
    adaptive_noise_variance = mean(mu_full.^2) / ...
        10^(adaptive_cfg.snr_dB / 10);

    initial_condition_trial = zeros(adaptive_cfg.num_monte_carlo, 1);
    extension_trial = zeros(adaptive_cfg.num_monte_carlo, 1);
    angle_error_trial = zeros(adaptive_cfg.num_monte_carlo, 3);
    nmse_trial = zeros(adaptive_cfg.num_monte_carlo, 3);

    for trial = 1:adaptive_cfg.num_monte_carlo
        rng(2000000 + trial, 'twister');
        z_initial = mu_initial + sqrt(adaptive_noise_variance) * ...
            randn(size(mu_initial));
        z_extra = mu_extra + sqrt(adaptive_noise_variance) * ...
            randn(size(mu_extra));
        Z_full = [reshape(z_initial, adaptive_cfg.num_pilots, ...
            adaptive_cfg.num_ris); reshape(z_extra, ...
            adaptive_cfg.num_pilots, adaptive_cfg.num_ris)];
        z_full = Z_full(:);

        [eta_initial_short, ~] = algorithm1_initialization( ...
            z_initial, initial_cycles, adaptive_cfg, ...
            coarse_dictionary, coarse_theta_deg, coarse_phi_deg);
        [~, J_initial_short] = training_model_and_jacobian( ...
            eta_initial_short, S_initial, B_initial, adaptive_cfg);
        singular_values = svd(J_initial_short, 'econ');
        initial_condition_trial(trial) = (max(singular_values) / ...
            max(min(singular_values), eps))^2;
        [eta_short, ~] = parametric_gn(z_initial, S_initial, ...
            B_initial, eta_initial_short, 1, adaptive_cfg);

        [eta_initial_full, ~] = algorithm1_initialization( ...
            z_full, full_cycles, full_cfg, coarse_dictionary, ...
            coarse_theta_deg, coarse_phi_deg);
        [eta_full, ~] = parametric_gn(z_full, S_full, B_full, ...
            eta_initial_full, 1, full_cfg);

        use_extension = initial_condition_trial(trial) > ...
            adaptive_jacobian_trigger || ...
            adaptive_pilot_condition(case_index) > adaptive_pilot_trigger;
        extension_trial(trial) = use_extension;
        if use_extension
            eta_adaptive = eta_full;
        else
            eta_adaptive = eta_short;
        end
        eta_candidates = {eta_short, eta_adaptive, eta_full};
        for method = 1:3
            eta_hat = eta_candidates{method};
            H_hat = eta_to_channel(eta_hat, cfg);
            [~, angle_error_trial(trial, method)] = matched_doa_error( ...
                eta_hat, cfg.theta_true_deg, cfg.phi_true_deg, cfg);
            nmse_trial(trial, method) = norm(H_hat - H_true, 'fro')^2 / ...
                max(norm(H_true, 'fro')^2, eps);
        end
    end

    adaptive_median_initial_condition(case_index) = median( ...
        initial_condition_trial);
    extension_probability(case_index) = mean(extension_trial);
    average_training_snapshots(case_index) = adaptive_cfg.num_pilots + ...
        adaptive_cfg.num_pilots * extension_probability(case_index);
    for method = 1:3
        adaptive_success(case_index, method) = mean( ...
            sqrt(angle_error_trial(:, method)) < ...
            adaptive_cfg.success_threshold_deg);
        adaptive_nmse_dB(case_index, method) = 10 * log10( ...
            max(mean(nmse_trial(:, method)), eps));
    end
    fprintf(['Adaptive rho %.3f, cond %.1f: extension %.3f, ', ...
        'success %.3f/%.3f/%.3f, NMSE %.2f/%.2f/%.2f dB\n'], ...
        rho_p, adaptive_pilot_condition(case_index), ...
        extension_probability(case_index), ...
        adaptive_success(case_index, :), adaptive_nmse_dB(case_index, :));
end

adaptive_table = table(adaptive_correlation(:), ...
    adaptive_pilot_condition, adaptive_median_initial_condition, ...
    extension_probability, average_training_snapshots, ...
    adaptive_success(:, 1), adaptive_success(:, 2), ...
    adaptive_success(:, 3), adaptive_nmse_dB(:, 1), ...
    adaptive_nmse_dB(:, 2), adaptive_nmse_dB(:, 3), ...
    'VariableNames', {'PilotCorrelation', 'InitialPilotGramCondition', ...
    'MedianInitialNormalCondition', 'ExtensionProbability', ...
    'AverageIntensitySnapshots', 'Fixed16_Success', ...
    'StabilityGuided_Success', 'Fixed32_Success', ...
    'Fixed16_Channel_NMSE_dB', 'StabilityGuided_Channel_NMSE_dB', ...
    'Fixed32_Channel_NMSE_dB'});
writetable(adaptive_table, fullfile(result_dir, ...
    'stability_guided_adaptive_training_summary.csv'));

fig = publication_figure([80, 80, 920, 350]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
adaptive_styles = {'^:', 'o-', 's--'};
adaptive_colors = {colors.gray, colors.blue, colors.orange};
nexttile;
for method = 1:3
    semilogx(adaptive_pilot_condition, adaptive_success(:, method), ...
        adaptive_styles{method}, 'Color', adaptive_colors{method});
    hold on;
end
grid on;
box on;
xlabel('Initial pilot Gram-matrix condition number');
ylabel('Successful recovery probability');
title('(a) Basin reliability');
legend(adaptive_method_names, 'Location', 'southwest');
ylim([0, 1.01]);
nexttile;
for method = 1:3
    semilogx(adaptive_pilot_condition, adaptive_nmse_dB(:, method), ...
        adaptive_styles{method}, 'Color', adaptive_colors{method});
    hold on;
end
grid on;
box on;
xlabel('Initial pilot Gram-matrix condition number');
ylabel('Channel NMSE (dB)');
title('(b) Estimation accuracy');
legend(adaptive_method_names, 'Location', 'northwest');
export_publication_figure(fig, result_dir, ...
    'fig_stability_guided_adaptive_training_accuracy');

fig = publication_figure([80, 80, 920, 350]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
loglog(adaptive_pilot_condition, adaptive_median_initial_condition, ...
    'd-', 'Color', colors.gray);
hold on;
yline(adaptive_jacobian_trigger, 'k--', '\kappa_{J,trigger}');
grid on;
box on;
xlabel('Initial pilot Gram-matrix condition number');
ylabel('Median \kappa(J_0^T J_0)');
title('(a) Stability trigger');
nexttile;
semilogx(adaptive_pilot_condition, average_training_snapshots, ...
    'o-', 'Color', colors.blue);
hold on;
yline(adaptive_cfg.num_pilots, ':', 'Fixed T_p=16', ...
    'Color', colors.gray);
yline(full_cfg.num_pilots, '--', 'Fixed T_p=32', ...
    'Color', colors.orange);
grid on;
box on;
xlabel('Initial pilot Gram-matrix condition number');
ylabel('Average intensity snapshots');
title('(b) Adaptive training overhead');
ylim([14, 34]);
export_publication_figure(fig, result_dir, ...
    'fig_stability_trigger_and_training_overhead');

save(fullfile(result_dir, ...
    'stability_guided_adaptive_training_results.mat'), ...
    'adaptive_cfg', 'full_cfg', 'adaptive_jacobian_trigger', ...
    'adaptive_pilot_trigger', ...
    'adaptive_correlation', 'adaptive_pilot_condition', ...
    'adaptive_median_initial_condition', 'extension_probability', ...
    'average_training_snapshots', 'adaptive_success', ...
    'adaptive_nmse_dB');

save(fullfile(result_dir, ...
    'stability_aware_practical_init_stress_results.mat'), ...
    'cfg', 'pilot_correlation', 'pilot_gram_condition', ...
    'median_initial_normal_condition', 'initial_doa_rmse_deg', ...
    'initial_channel_nmse_dB', 'success_probability', ...
    'catastrophic_failure_probability', 'doa_rmse_deg', ...
    'channel_nmse_dB', 'average_iterations', 'average_rejections', ...
    'stability_activation_probability');

fprintf('\nAll outputs saved to:\n%s\n', result_dir);

%% Local functions
function [S, B, pilot_cycles] = make_four_phase_training(cfg, rho_p, seed)
    old_rng = rng;
    rng(seed, 'twister');
    qpsk = [1 + 1j; -1 + 1j; -1 - 1j; 1 - 1j] / sqrt(2);
    first_pilot = qpsk(randi(4, cfg.num_cycles, 1));
    orthogonal_seed = qpsk(randi(4, cfg.num_cycles, 1));
    orthogonal_pilot = orthogonal_seed - first_pilot * ...
        ((first_pilot' * orthogonal_seed) / (first_pilot' * first_pilot));
    if norm(orthogonal_pilot) < 1e-8
        orthogonal_seed = exp(1j * 2 * pi * ...
            (0:cfg.num_cycles - 1).' / cfg.num_cycles);
        orthogonal_pilot = orthogonal_seed - first_pilot * ...
            ((first_pilot' * orthogonal_seed) / ...
            (first_pilot' * first_pilot));
    end
    orthogonal_pilot = sqrt(cfg.num_cycles) * orthogonal_pilot / ...
        norm(orthogonal_pilot);
    second_pilot = rho_p * first_pilot + sqrt(1 - rho_p^2) * ...
        orthogonal_pilot;
    pilot_cycles = [first_pilot, second_pilot];
    assert(rank(pilot_cycles) == cfg.num_users, ...
        'All pilot matrices in this experiment must remain full rank.');

    S = repelem(pilot_cycles, cfg.num_reference_states, 1);
    reference_phase = repmat((0:cfg.num_reference_states - 1).' * ...
        pi / 2, cfg.num_cycles, 1);
    space_phase = 2 * pi * (0:cfg.num_ris - 1) / cfg.num_ris;
    B = cfg.reference_amplitude * exp( ...
        1j * (reference_phase + space_phase));
    rng(old_rng);
end

function [initial_cycles, complementary_cycles] = ...
        make_complementary_pilot_cycles(cfg, rho_p, seed)
    old_rng = rng;
    rng(seed, 'twister');
    qpsk = [1 + 1j; -1 + 1j; -1 - 1j; 1 - 1j] / sqrt(2);
    first_pilot = qpsk(randi(4, cfg.num_cycles, 1));
    orthogonal_seed = qpsk(randi(4, cfg.num_cycles, 1));
    orthogonal_pilot = orthogonal_seed - first_pilot * ...
        ((first_pilot' * orthogonal_seed) / (first_pilot' * first_pilot));
    if norm(orthogonal_pilot) < 1e-8
        orthogonal_seed = exp(1j * 2 * pi * ...
            (0:cfg.num_cycles - 1).' / cfg.num_cycles);
        orthogonal_pilot = orthogonal_seed - first_pilot * ...
            ((first_pilot' * orthogonal_seed) / ...
            (first_pilot' * first_pilot));
    end
    orthogonal_pilot = sqrt(cfg.num_cycles) * orthogonal_pilot / ...
        norm(orthogonal_pilot);
    second_initial = rho_p * first_pilot + sqrt(1 - rho_p^2) * ...
        orthogonal_pilot;
    second_complementary = -rho_p * first_pilot + ...
        sqrt(1 - rho_p^2) * orthogonal_pilot;
    initial_cycles = [first_pilot, second_initial];
    complementary_cycles = [first_pilot, second_complementary];
    assert(rank(initial_cycles) == cfg.num_users, ...
        'The initial pilot matrix must remain full rank.');
    assert(cond([initial_cycles; complementary_cycles]' * ...
        [initial_cycles; complementary_cycles]) < 1 + 1e-10, ...
        'The complementary block should orthogonalize the combined pilots.');
    rng(old_rng);
end

function [S, B] = build_four_phase_training(pilot_cycles, cfg)
    num_cycles = size(pilot_cycles, 1);
    S = repelem(pilot_cycles, cfg.num_reference_states, 1);
    reference_phase = repmat((0:cfg.num_reference_states - 1).' * ...
        pi / 2, num_cycles, 1);
    space_phase = 2 * pi * (0:cfg.num_ris - 1) / cfg.num_ris;
    B = cfg.reference_amplitude * exp( ...
        1j * (reference_phase + space_phase));
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
    Z = reshape(z, cfg.num_pilots, cfg.num_ris);
    space_phase = 2 * pi * (0:cfg.num_ris - 1) / cfg.num_ris;
    field_cycles = zeros(cfg.num_cycles, cfg.num_ris);
    for cycle = 1:cfg.num_cycles
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

function [mu, J] = training_model_and_jacobian(eta, S, B, cfg)
    L = cfg.num_users;
    theta_deg = eta(1:L) * 180 / pi;
    phi_deg = eta(L + 1:2 * L) * 180 / pi;
    alpha = eta(2 * L + 1:3 * L) + ...
        1j * eta(3 * L + 1:4 * L);
    field_derivative = zeros(numel(B), 4 * L);
    field = zeros(size(B));
    for user = 1:L
        [v, dv_theta, dv_phi] = steering_vector_with_derivatives( ...
            theta_deg(user), phi_deg(user), cfg);
        field = field + S(:, user) * (alpha(user) * v).';
        derivatives = {alpha(user) * dv_theta, ...
            alpha(user) * dv_phi, v, 1j * v};
        columns = [user, L + user, 2 * L + user, 3 * L + user];
        for derivative_index = 1:4
            derivative_matrix = S(:, user) * ...
                derivatives{derivative_index}.';
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
    rejections = 0;
    stability_activations = 0;
    for iteration = 1:cfg.gn_max_iter
        singular_values = svd(J, 'econ');
        sigma_max = max(singular_values);
        sigma_min = min(singular_values);
        if method == 1
            lambda = lambda_lm;
        else
            lambda_stability = max((sigma_max^2 - ...
                cfg.gn_target_condition * sigma_min^2) / ...
                (cfg.gn_target_condition - 1), 0);
            if lambda_stability > lambda_lm
                stability_activations = stability_activations + 1;
            end
            lambda = max(lambda_lm, lambda_stability);
        end
        increment = (J.' * J + (lambda + 1e-10) * ...
            eye(size(J, 2))) \ (J.' * residual);
        eta_trial = project_eta(eta + increment, cfg);
        [mu_trial, J_trial] = training_model_and_jacobian( ...
            eta_trial, S, B, cfg);
        residual_trial = z - mu_trial;
        cost_trial = 0.5 * norm(residual_trial)^2;
        if cost_trial < cost
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
        if norm(increment) <= cfg.gn_tolerance * (norm(eta) + 1)
            break;
        end
    end
    information.iterations = iteration;
    information.rejections = rejections;
    information.stability_activations = stability_activations;
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

function fig = publication_figure(position)
    fig = figure('Color', 'w', 'Position', position);
    set(groot, 'defaultAxesFontName', 'Times New Roman');
    set(groot, 'defaultTextFontName', 'Times New Roman');
    set(groot, 'defaultAxesFontSize', 11);
    set(groot, 'defaultLineLineWidth', 1.5);
    set(groot, 'defaultLineMarkerSize', 7);
end

function export_publication_figure(fig, result_dir, base_name)
    set(findall(fig, '-property', 'LineWidth'), 'LineWidth', 1.5);
    exportgraphics(fig, fullfile(result_dir, [base_name, '.eps']), ...
        'ContentType', 'vector');
    exportgraphics(fig, fullfile(result_dir, [base_name, '.png']), ...
        'Resolution', 300);
    savefig(fig, fullfile(result_dir, [base_name, '.fig']));
end
