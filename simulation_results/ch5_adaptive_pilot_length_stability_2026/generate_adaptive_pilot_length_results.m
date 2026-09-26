%% Stability-guided adaptive pilot-length selection
clc;
clear;
close all;

script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

rng(20260925, 'twister');
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
cfg.min_cycles = 2;
cfg.max_cycles = 10;
cfg.path_gains = [exp(1j * 0.35); 0.85 * exp(-1j * 0.55)];
cfg.theta_true_deg = [19; 21];
cfg.phi_true_deg = [5; 5];
cfg.target_doa_crlb_deg = 1.0;
cfg.target_channel_crlb_nmse_dB = -25;
cfg.target_normal_condition = 5e3;
cfg.success_threshold_deg = 1;
cfg.gn_max_iter = 40;
cfg.gn_tolerance = 1e-6;
cfg.gn_initial_damping = 1e-3;
cfg.coarse_theta_grid_deg = 0:0.5:60;
cfg.coarse_phi_grid_deg = -30:1:30;
cfg.num_monte_carlo = 300;
if quick_check
    cfg.num_monte_carlo = 20;
end

snr_dB = 0:2:20;
num_snr = numel(snr_dB);
method_names = {'Fixed short training', ...
    'Stability-guided adaptive', 'Fixed long training'};
num_methods = numel(method_names);

colors.gray = [0.35, 0.35, 0.35];
colors.blue = [0.0000, 0.4470, 0.7410];
colors.orange = [0.8500, 0.3250, 0.0980];

[coarse_dictionary, coarse_theta_deg, coarse_phi_deg] = ...
    make_coarse_dictionary(cfg);
eta_true = pack_eta(cfg.theta_true_deg, cfg.phi_true_deg, cfg.path_gains);
H_true = eta_to_channel(eta_true, cfg);
pilot_cycles_full = make_nested_qpsk_pilots(cfg);

cfg_full = cfg;
cfg_full.num_cycles = cfg.max_cycles;
cfg_full.num_pilots = cfg.max_cycles * cfg.num_reference_states;
[S_full, B_full] = build_four_phase_training( ...
    pilot_cycles_full, cfg_full);
[mu_full, ~] = training_model_and_jacobian( ...
    eta_true, S_full, B_full, cfg_full);

channel_nmse_dB = zeros(num_snr, num_methods);
doa_rmse_deg = zeros(num_snr, num_methods);
success_probability = zeros(num_snr, num_methods);
average_selected_pilots = zeros(num_snr, 1);
median_selected_pilots = zeros(num_snr, 1);
max_length_probability = zeros(num_snr, 1);
average_predicted_crlb_deg = zeros(num_snr, 1);
average_predicted_channel_nmse_dB = zeros(num_snr, 1);
selection_histogram = zeros(num_snr, cfg.max_cycles - cfg.min_cycles + 1);

fprintf('\nStability-guided adaptive pilot-length experiment\n');
fprintf(['RIS %d x %d, users %d, cycle range %d--%d, ', ...
    'trials %d\n'], cfg.ris_rows, cfg.ris_cols, cfg.num_users, ...
    cfg.min_cycles, cfg.max_cycles, cfg.num_monte_carlo);

for snr_index = 1:num_snr
    noise_variance = mean(mu_full.^2) / 10^(snr_dB(snr_index) / 10);
    nmse_trial = zeros(cfg.num_monte_carlo, num_methods);
    angle_error_trial = zeros(cfg.num_monte_carlo, num_methods);
    selected_pilots_trial = zeros(cfg.num_monte_carlo, 1);
    selected_crlb_trial = zeros(cfg.num_monte_carlo, 1);
    selected_channel_bound_trial = zeros(cfg.num_monte_carlo, 1);

    for trial = 1:cfg.num_monte_carlo
        rng(3000000 + 10000 * snr_index + trial, 'twister');
        z_full = mu_full + sqrt(noise_variance) * randn(size(mu_full));
        Z_full = reshape(z_full, cfg_full.num_pilots, cfg.num_ris);

        eta_initial_store = cell(cfg.max_cycles, 1);
        z_store = cell(cfg.max_cycles, 1);
        S_store = cell(cfg.max_cycles, 1);
        B_store = cell(cfg.max_cycles, 1);
        cfg_store = cell(cfg.max_cycles, 1);
        predicted_crlb_store = inf(cfg.max_cycles, 1);
        predicted_channel_bound_store = inf(cfg.max_cycles, 1);
        condition_store = inf(cfg.max_cycles, 1);
        selected_cycles = cfg.max_cycles;

        for num_cycles = cfg.min_cycles:cfg.max_cycles
            current_cfg = cfg;
            current_cfg.num_cycles = num_cycles;
            current_cfg.num_pilots = num_cycles * ...
                cfg.num_reference_states;
            pilot_cycles = pilot_cycles_full(1:num_cycles, :);
            [S_current, B_current] = build_four_phase_training( ...
                pilot_cycles, current_cfg);
            z_current = Z_full(1:current_cfg.num_pilots, :);
            z_current = z_current(:);
            [eta_initial, ~] = algorithm1_initialization( ...
                z_current, pilot_cycles, current_cfg, ...
                coarse_dictionary, coarse_theta_deg, coarse_phi_deg);
            [~, J_initial] = training_model_and_jacobian( ...
                eta_initial, S_current, B_current, current_cfg);
            [predicted_crlb, predicted_channel_bound, ...
                normal_condition] = fim_accuracy_metrics( ...
                eta_initial, J_initial, noise_variance, current_cfg);

            eta_initial_store{num_cycles} = eta_initial;
            z_store{num_cycles} = z_current;
            S_store{num_cycles} = S_current;
            B_store{num_cycles} = B_current;
            cfg_store{num_cycles} = current_cfg;
            predicted_crlb_store(num_cycles) = predicted_crlb;
            predicted_channel_bound_store(num_cycles) = ...
                predicted_channel_bound;
            condition_store(num_cycles) = normal_condition;

            if predicted_crlb <= cfg.target_doa_crlb_deg && ...
                    predicted_channel_bound <= ...
                    cfg.target_channel_crlb_nmse_dB && ...
                    normal_condition <= cfg.target_normal_condition
                selected_cycles = num_cycles;
                break;
            end
        end

        short_cycles = cfg.min_cycles;
        long_cycles = cfg.max_cycles;
        required_cycles = unique([short_cycles, selected_cycles, long_cycles]);
        eta_hat_store = cell(cfg.max_cycles, 1);
        for cycle_index = 1:numel(required_cycles)
            num_cycles = required_cycles(cycle_index);
            if isempty(eta_initial_store{num_cycles})
                current_cfg = cfg;
                current_cfg.num_cycles = num_cycles;
                current_cfg.num_pilots = num_cycles * ...
                    cfg.num_reference_states;
                pilot_cycles = pilot_cycles_full(1:num_cycles, :);
                [S_current, B_current] = build_four_phase_training( ...
                    pilot_cycles, current_cfg);
                z_current = Z_full(1:current_cfg.num_pilots, :);
                z_current = z_current(:);
                [eta_initial_store{num_cycles}, ~] = ...
                    algorithm1_initialization(z_current, pilot_cycles, ...
                    current_cfg, coarse_dictionary, coarse_theta_deg, ...
                    coarse_phi_deg);
                z_store{num_cycles} = z_current;
                S_store{num_cycles} = S_current;
                B_store{num_cycles} = B_current;
                cfg_store{num_cycles} = current_cfg;
            end
            [eta_hat_store{num_cycles}, ~] = parametric_lm( ...
                z_store{num_cycles}, S_store{num_cycles}, ...
                B_store{num_cycles}, eta_initial_store{num_cycles}, ...
                cfg_store{num_cycles});
        end

        eta_candidates = {eta_hat_store{short_cycles}, ...
            eta_hat_store{selected_cycles}, eta_hat_store{long_cycles}};
        for method = 1:num_methods
            eta_hat = eta_candidates{method};
            H_hat = eta_to_channel(eta_hat, cfg);
            [~, angle_error_trial(trial, method)] = matched_doa_error( ...
                eta_hat, cfg.theta_true_deg, cfg.phi_true_deg, cfg);
            nmse_trial(trial, method) = norm(H_hat - H_true, 'fro')^2 / ...
                max(norm(H_true, 'fro')^2, eps);
        end
        selected_pilots_trial(trial) = selected_cycles * ...
            cfg.num_reference_states;
        selected_crlb_trial(trial) = ...
            predicted_crlb_store(selected_cycles);
        selected_channel_bound_trial(trial) = ...
            predicted_channel_bound_store(selected_cycles);
        selection_histogram(snr_index, ...
            selected_cycles - cfg.min_cycles + 1) = ...
            selection_histogram(snr_index, ...
            selected_cycles - cfg.min_cycles + 1) + 1;
    end

    for method = 1:num_methods
        channel_nmse_dB(snr_index, method) = 10 * log10( ...
            max(mean(nmse_trial(:, method)), eps));
        doa_rmse_deg(snr_index, method) = sqrt(mean( ...
            angle_error_trial(:, method)));
        success_probability(snr_index, method) = mean( ...
            sqrt(angle_error_trial(:, method)) < ...
            cfg.success_threshold_deg);
    end
    average_selected_pilots(snr_index) = mean(selected_pilots_trial);
    median_selected_pilots(snr_index) = median(selected_pilots_trial);
    max_length_probability(snr_index) = mean( ...
        selected_pilots_trial == cfg_full.num_pilots);
    average_predicted_crlb_deg(snr_index) = mean(selected_crlb_trial);
    average_predicted_channel_nmse_dB(snr_index) = mean( ...
        selected_channel_bound_trial);
    selection_histogram(snr_index, :) = ...
        selection_histogram(snr_index, :) / cfg.num_monte_carlo;

    fprintf(['SNR %4.1f dB: selected T_p %.2f, NMSE ', ...
        '[%.2f %.2f %.2f] dB, success [%.3f %.3f %.3f]\n'], ...
        snr_dB(snr_index), average_selected_pilots(snr_index), ...
        channel_nmse_dB(snr_index, :), ...
        success_probability(snr_index, :));
end

%% Save results
summary_table = table(snr_dB(:), channel_nmse_dB(:, 1), ...
    channel_nmse_dB(:, 2), channel_nmse_dB(:, 3), ...
    doa_rmse_deg(:, 1), doa_rmse_deg(:, 2), doa_rmse_deg(:, 3), ...
    success_probability(:, 1), success_probability(:, 2), ...
    success_probability(:, 3), average_selected_pilots, ...
    median_selected_pilots, max_length_probability, ...
    average_predicted_crlb_deg, average_predicted_channel_nmse_dB, ...
    'VariableNames', {'SNR_dB', 'FixedShort_Channel_NMSE_dB', ...
    'Adaptive_Channel_NMSE_dB', 'FixedLong_Channel_NMSE_dB', ...
    'FixedShort_DOA_RMSE_deg', 'Adaptive_DOA_RMSE_deg', ...
    'FixedLong_DOA_RMSE_deg', 'FixedShort_Success', ...
    'Adaptive_Success', 'FixedLong_Success', ...
    'AverageSelectedPilots', 'MedianSelectedPilots', ...
    'MaxLengthProbability', 'AveragePredictedDOA_CRLB_deg', ...
    'AveragePredictedChannel_CRLB_NMSE_dB'});
writetable(summary_table, fullfile(result_dir, ...
    'adaptive_pilot_length_summary.csv'));

histogram_table = array2table(selection_histogram, ...
    'VariableNames', compose('Tp_%d', ...
    (cfg.min_cycles:cfg.max_cycles) * cfg.num_reference_states));
histogram_table = addvars(histogram_table, snr_dB(:), ...
    'Before', 1, 'NewVariableNames', 'SNR_dB');
writetable(histogram_table, fullfile(result_dir, ...
    'adaptive_pilot_length_distribution.csv'));

%% Requested two-panel figure
fig = publication_figure([80, 80, 920, 350]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
styles = {'^:', 'o-', 's--'};
method_colors = {colors.gray, colors.blue, colors.orange};
nexttile;
for method = 1:num_methods
    plot(snr_dB, channel_nmse_dB(:, method), styles{method}, ...
        'Color', method_colors{method});
    hold on;
end
grid on;
box on;
xlabel('Training SNR (dB)');
ylabel('Channel NMSE (dB)');
title('(a) Estimation accuracy');
legend(method_names, 'Location', 'southwest');
nexttile;
plot(snr_dB, average_selected_pilots, 'o-', 'Color', colors.blue);
hold on;
yline(cfg.min_cycles * cfg.num_reference_states, ':', ...
    'Fixed short training', 'Color', colors.gray);
yline(cfg.max_cycles * cfg.num_reference_states, '--', ...
    'Fixed long training', 'Color', colors.orange);
grid on;
box on;
xlabel('Training SNR (dB)');
ylabel('Average selected pilot measurements');
title('(b) Adaptive training overhead');
ylim([cfg.min_cycles * cfg.num_reference_states - 2, ...
    cfg.max_cycles * cfg.num_reference_states + 2]);
export_publication_figure(fig, result_dir, ...
    'fig_adaptive_pilot_length_nmse_overhead');

save(fullfile(result_dir, 'adaptive_pilot_length_results.mat'), ...
    'cfg', 'snr_dB', 'channel_nmse_dB', 'doa_rmse_deg', ...
    'success_probability', 'average_selected_pilots', ...
    'median_selected_pilots', 'max_length_probability', ...
    'average_predicted_crlb_deg', ...
    'average_predicted_channel_nmse_dB', 'selection_histogram');

fprintf('\nAll outputs saved to:\n%s\n', result_dir);

%% Local functions
function pilot_cycles = make_nested_qpsk_pilots(cfg)
    pattern = repmat([1; -1; 1j; -1j], ...
        ceil(cfg.max_cycles / 4), 1);
    pattern = pattern(1:cfg.max_cycles);
    first_pilot = ones(cfg.max_cycles, 1);
    second_pilot = first_pilot .* pattern;
    pilot_cycles = [first_pilot, second_pilot];
    for num_cycles = cfg.min_cycles:cfg.max_cycles
        assert(rank(pilot_cycles(1:num_cycles, :)) == cfg.num_users, ...
            'Every candidate pilot prefix must have full column rank.');
    end
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

function [crlb_rmse_deg, channel_nmse_bound_dB, normal_condition] = ...
        fim_accuracy_metrics(eta, J, noise_variance, cfg)
    L = cfg.num_users;
    F = (J' * J) / max(noise_variance, eps);
    singular_values = svd(J, 'econ');
    normal_condition = (max(singular_values) / ...
        max(min(singular_values), eps))^2;
    F_pp = F(1:2 * L, 1:2 * L);
    F_pn = F(1:2 * L, 2 * L + 1:end);
    F_nn = F(2 * L + 1:end, 2 * L + 1:end);
    F_doa = F_pp - F_pn * pinv(F_nn) * F_pn';
    F_doa = 0.5 * (F_doa + F_doa');
    covariance_bound = pinv(F_doa);
    bound_radians = sqrt(max(real(trace(covariance_bound)) / ...
        (2 * L), 0));
    crlb_rmse_deg = bound_radians * 180 / pi;

    parameter_covariance = pinv(F);
    channel_jacobian = channel_parameter_jacobian(eta, cfg);
    predicted_channel_mse = real(trace(channel_jacobian * ...
        parameter_covariance * channel_jacobian'));
    H_estimate = eta_to_channel(eta, cfg);
    channel_nmse_bound_dB = 10 * log10(max( ...
        predicted_channel_mse / max(norm(H_estimate, 'fro')^2, eps), eps));
end

function D = channel_parameter_jacobian(eta, cfg)
    L = cfg.num_users;
    theta_deg = eta(1:L) * 180 / pi;
    phi_deg = eta(L + 1:2 * L) * 180 / pi;
    alpha = eta(2 * L + 1:3 * L) + ...
        1j * eta(3 * L + 1:4 * L);
    D = zeros(cfg.num_ris * L, 4 * L);
    for user = 1:L
        [v, dv_theta, dv_phi] = steering_vector_with_derivatives( ...
            theta_deg(user), phi_deg(user), cfg);
        rows = (user - 1) * cfg.num_ris + (1:cfg.num_ris);
        D(rows, user) = alpha(user) * dv_theta;
        D(rows, L + user) = alpha(user) * dv_phi;
        D(rows, 2 * L + user) = v;
        D(rows, 3 * L + user) = 1j * v;
    end
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

function [eta, information] = parametric_lm(z, S, B, eta_initial, cfg)
    eta = eta_initial;
    lambda_lm = cfg.gn_initial_damping;
    [mu, J] = training_model_and_jacobian(eta, S, B, cfg);
    residual = z - mu;
    cost = 0.5 * norm(residual)^2;
    rejections = 0;
    for iteration = 1:cfg.gn_max_iter
        increment = (J.' * J + (lambda_lm + 1e-10) * ...
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
