%% Stability-guided adaptive training with a fixed spatial four-phase code
clc;
clear;
close all;

script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

rng(20260927, 'twister');
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');
quick_check = strcmp(getenv('HOLO_QUICK_CHECK'), '1');
use_parallel = license('test', 'Distrib_Computing_Toolbox');
max_parallel_workers = 0;
if use_parallel
    max_parallel_workers = 4;
end

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
cfg.residual_discrepancy_threshold = 1.10;
cfg.success_threshold_deg = 1;
cfg.gn_max_iter = 40;
cfg.gn_tolerance = 1e-6;
cfg.gn_initial_damping = 1e-3;
cfg.coarse_theta_grid_deg = 0:0.5:60;
cfg.coarse_phi_grid_deg = -30:1:30;
cfg.num_monte_carlo = 300;
cfg.data_trials_per_channel = 200;
if quick_check
    cfg.num_monte_carlo = 20;
end

snr_dB = 0:2:20;
num_snr = numel(snr_dB);
method_names = {'Fixed short', 'Residual-adaptive', ...
    'Stability-guided', 'Oracle stopping', 'Fixed long'};
num_methods = numel(method_names);
short_method = 1;
residual_method = 2;
proposed_method = 3;
oracle_method = 4;
long_method = 5;

colors.gray = [0.35, 0.35, 0.35];
colors.purple = [0.4940, 0.1840, 0.5560];
colors.blue = [0.0000, 0.4470, 0.7410];
colors.black = [0.10, 0.10, 0.10];
colors.orange = [0.8500, 0.3250, 0.0980];

[coarse_dictionary, coarse_theta_deg, coarse_phi_deg] = ...
    make_coarse_dictionary(cfg);
eta_true = pack_eta(cfg.theta_true_deg, cfg.phi_true_deg, cfg.path_gains);
H_true = eta_to_channel(eta_true, cfg);
pilot_cycles_full = make_nested_qpsk_pilots(cfg);

cfg_full = cfg;
cfg_full.num_cycles = cfg.max_cycles;
cfg_full.num_pilots = cfg.max_cycles * cfg.num_reference_states;
[S_full, B_full] = build_spatial_reference_training(pilot_cycles_full, cfg_full);
[mu_full, ~] = training_model_and_jacobian( ...
    eta_true, S_full, B_full, cfg_full);

channel_nmse_dB = zeros(num_snr, num_methods);
channel_nmse_ci95_dB = zeros(num_snr, num_methods);
doa_rmse_deg = zeros(num_snr, num_methods);
success_probability = zeros(num_snr, num_methods);
average_selected_pilots = zeros(num_snr, num_methods);
median_selected_pilots = zeros(num_snr, num_methods);
selection_histogram = zeros(num_snr, num_methods, ...
    cfg.max_cycles - cfg.min_cycles + 1);
average_predicted_channel_nmse_dB = zeros(num_snr, 1);
end_to_end_ser = zeros(num_snr, 4);
qpsk = exp(1j * (pi / 4 + (0:3) * pi / 2)).';
data_symbols = symbol_codebook(qpsk, cfg.num_users);
data_reference = cfg.reference_amplitude * exp(1j * ...
    (pi / 2) * mod((0:cfg.num_ris - 1).', 4));
data_mu_true = abs(H_true * data_symbols + data_reference).^2;
data_snr_dB = 0;
data_noise_variance = mean(data_mu_true(:).^2) / 10^(data_snr_dB / 10);

fprintf('\nPractical adaptive-training baseline comparison\n');
fprintf(['Surface %d x %d, users %d, pilot symbols %d--%d, ', ...
    'trials %d\n'], cfg.ris_rows, cfg.ris_cols, cfg.num_users, ...
    cfg.min_cycles * cfg.num_reference_states, ...
    cfg.max_cycles * cfg.num_reference_states, cfg.num_monte_carlo);
if use_parallel && isempty(gcp('nocreate'))
    parpool('local', max_parallel_workers);
end

for snr_index = 1:num_snr
    noise_variance = mean(mu_full.^2) / 10^(snr_dB(snr_index) / 10);
    nmse_trial = zeros(cfg.num_monte_carlo, num_methods);
    angle_error_trial = zeros(cfg.num_monte_carlo, num_methods);
    selected_cycles_trial = zeros(cfg.num_monte_carlo, num_methods);
    selected_bound_trial = zeros(cfg.num_monte_carlo, 1);
    ser_trial = zeros(cfg.num_monte_carlo, 4);

    parfor (trial = 1:cfg.num_monte_carlo, max_parallel_workers)
        rng(4000000 + 10000 * snr_index + trial, 'twister');
        z_full = mu_full + sqrt(noise_variance) * randn(size(mu_full));
        Z_full = reshape(z_full, cfg_full.num_pilots, cfg.num_ris);

        eta_hat_store = cell(cfg.max_cycles, 1);
        predicted_bound_store = inf(cfg.max_cycles, 1);
        selected_residual = cfg.max_cycles;
        selected_proposed = cfg.max_cycles;
        selected_oracle = cfg.max_cycles;
        residual_found = false;
        proposed_found = false;
        oracle_found = false;

        for num_cycles = cfg.min_cycles:cfg.max_cycles
            current_cfg = cfg;
            current_cfg.num_cycles = num_cycles;
            current_cfg.num_pilots = num_cycles * ...
                cfg.num_reference_states;
            current_cfg.current_noise_variance = noise_variance;
            pilot_cycles = pilot_cycles_full(1:current_cfg.num_pilots, :);
            [S_current, B_current] = build_spatial_reference_training( ...
                pilot_cycles, current_cfg);
            z_current = Z_full(1:current_cfg.num_pilots, :);
            z_current = z_current(:);

            [eta_initial, ~] = algorithm1_initialization( ...
                z_current, pilot_cycles, current_cfg, ...
                coarse_dictionary, coarse_theta_deg, coarse_phi_deg);
            [~, J_initial] = training_model_and_jacobian( ...
                eta_initial, S_current, B_current, current_cfg);
            [predicted_doa_crlb, predicted_channel_bound, ...
                normal_condition] = fim_accuracy_metrics( ...
                eta_initial, J_initial, noise_variance, current_cfg);
            predicted_bound_store(num_cycles) = predicted_channel_bound;

            [eta_hat, solver_info] = parametric_lm( ...
                z_current, S_current, B_current, eta_initial, current_cfg);
            eta_hat_store{num_cycles} = eta_hat;

            if ~residual_found && ...
                    solver_info.reduced_residual <= ...
                    cfg.residual_discrepancy_threshold
                selected_residual = num_cycles;
                residual_found = true;
            end

            if ~proposed_found && ...
                    predicted_doa_crlb <= cfg.target_doa_crlb_deg && ...
                    predicted_channel_bound <= ...
                    cfg.target_channel_crlb_nmse_dB && ...
                    normal_condition <= cfg.target_normal_condition
                selected_proposed = num_cycles;
                proposed_found = true;
            end

            H_hat_current = eta_to_channel(eta_hat, cfg);
            current_nmse_dB = 10 * log10(max( ...
                norm(H_hat_current - H_true, 'fro')^2 / ...
                max(norm(H_true, 'fro')^2, eps), eps));
            if ~oracle_found && current_nmse_dB <= ...
                    cfg.target_channel_crlb_nmse_dB
                selected_oracle = num_cycles;
                oracle_found = true;
            end
        end

        selected_cycles = [cfg.min_cycles, selected_residual, ...
            selected_proposed, selected_oracle, cfg.max_cycles];
        for method = 1:num_methods
            eta_hat = eta_hat_store{selected_cycles(method)};
            H_hat = eta_to_channel(eta_hat, cfg);
            nmse_trial(trial, method) = norm(H_hat - H_true, 'fro')^2 / ...
                max(norm(H_true, 'fro')^2, eps);
            [~, angle_error_trial(trial, method)] = matched_doa_error( ...
                eta_hat, cfg.theta_true_deg, cfg.phi_true_deg, cfg);
            selected_cycles_trial(trial, method) = selected_cycles(method);
        end
        compared_methods = [short_method, proposed_method, long_method];
        error_count = zeros(1, 4);
        for data_trial = 1:cfg.data_trials_per_channel
            symbol_index = randi(size(data_symbols, 2));
            z_data = data_mu_true(:, symbol_index) + ...
                sqrt(data_noise_variance) * randn(cfg.num_ris, 1);
            for compared = 1:3
                eta_data = eta_hat_store{selected_cycles(compared_methods(compared))};
                H_data = eta_to_channel(eta_data, cfg);
                mu_data = abs(H_data * data_symbols + data_reference).^2;
                [~, detected] = min(sum((mu_data - z_data).^2, 1));
                error_count(compared) = error_count(compared) + ...
                    sum(data_symbols(:, detected) ~= data_symbols(:, symbol_index));
            end
            [~, detected] = min(sum((data_mu_true - z_data).^2, 1));
            error_count(4) = error_count(4) + ...
                sum(data_symbols(:, detected) ~= data_symbols(:, symbol_index));
        end
        ser_trial(trial, :) = error_count / ...
            (cfg.data_trials_per_channel * cfg.num_users);
        selected_bound_trial(trial) = ...
            predicted_bound_store(selected_proposed);
    end

    for method = 1:num_methods
        mean_nmse = mean(nmse_trial(:, method));
        se_nmse = std(nmse_trial(:, method)) / ...
            sqrt(cfg.num_monte_carlo);
        channel_nmse_dB(snr_index, method) = ...
            10 * log10(max(mean_nmse, eps));
        upper_nmse = max(mean_nmse + 1.96 * se_nmse, eps);
        channel_nmse_ci95_dB(snr_index, method) = ...
            10 * log10(upper_nmse) - channel_nmse_dB(snr_index, method);
        doa_rmse_deg(snr_index, method) = sqrt(mean( ...
            angle_error_trial(:, method)));
        success_probability(snr_index, method) = mean( ...
            sqrt(angle_error_trial(:, method)) < ...
            cfg.success_threshold_deg);
        average_selected_pilots(snr_index, method) = mean( ...
            selected_cycles_trial(:, method)) * cfg.num_reference_states;
        median_selected_pilots(snr_index, method) = median( ...
            selected_cycles_trial(:, method)) * cfg.num_reference_states;
        for num_cycles = cfg.min_cycles:cfg.max_cycles
            selection_histogram(snr_index, method, ...
                num_cycles - cfg.min_cycles + 1) = mean( ...
                selected_cycles_trial(:, method) == num_cycles);
        end
    end
    average_predicted_channel_nmse_dB(snr_index) = mean( ...
        selected_bound_trial);
    end_to_end_ser(snr_index, :) = mean(ser_trial, 1);

    fprintf(['SNR %4.1f dB: NMSE [', repmat('%7.2f ', 1, num_methods), ...
        '] dB, average T_p [', repmat('%5.1f ', 1, num_methods), ']\n'], ...
        snr_dB(snr_index), channel_nmse_dB(snr_index, :), ...
        average_selected_pilots(snr_index, :));
end

%% Save numerical results
summary_table = table(snr_dB(:), 'VariableNames', {'SNR_dB'});
for method = 1:num_methods
    safe_name = matlab.lang.makeValidName(method_names{method});
    summary_table.([safe_name, '_Channel_NMSE_dB']) = ...
        channel_nmse_dB(:, method);
    summary_table.([safe_name, '_NMSE_CI95_dB']) = ...
        channel_nmse_ci95_dB(:, method);
    summary_table.([safe_name, '_DOA_RMSE_deg']) = ...
        doa_rmse_deg(:, method);
    summary_table.([safe_name, '_Success']) = ...
        success_probability(:, method);
    summary_table.([safe_name, '_Average_Tp']) = ...
        average_selected_pilots(:, method);
    summary_table.([safe_name, '_Median_Tp']) = ...
        median_selected_pilots(:, method);
end
summary_table.ProposedPredictedChannelCRLBNMSE_dB = ...
    average_predicted_channel_nmse_dB;
writetable(summary_table, fullfile(result_dir, ...
    'stability_guided_training_baseline_summary.csv'));

histogram_rows = [];
for snr_index = 1:num_snr
    for method = 1:num_methods
        for num_cycles = cfg.min_cycles:cfg.max_cycles
            histogram_rows = [histogram_rows; ...
                snr_dB(snr_index), method, ...
                num_cycles * cfg.num_reference_states, ...
                selection_histogram(snr_index, method, ...
                num_cycles - cfg.min_cycles + 1)]; %#ok<AGROW>
        end
    end
end
histogram_table = array2table(histogram_rows, 'VariableNames', ...
    {'SNR_dB', 'MethodIndex', 'SelectedTp', 'Probability'});
histogram_table.Method = categorical(method_names(histogram_table.MethodIndex)).';
histogram_table = movevars(histogram_table, 'Method', 'After', 'MethodIndex');
writetable(histogram_table, fullfile(result_dir, ...
    'stability_guided_training_selection_distribution.csv'));

%% Two-panel publication figure
fig = publication_figure([80, 80, 940, 350]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
styles = {'^:', 'd-.', 'o-', 'x--', 's--'};
method_colors = {colors.gray, colors.purple, colors.blue, ...
    colors.black, colors.orange};

nexttile;
for method = 1:num_methods
    plot(snr_dB, channel_nmse_dB(:, method), styles{method}, ...
        'Color', method_colors{method});
    hold on;
end
yline(cfg.target_channel_crlb_nmse_dB, 'k:', ...
    'Target', 'HandleVisibility', 'off');
grid on;
box on;
xlabel('Training SNR (dB)');
ylabel('Channel NMSE (dB)');
title('(a) Estimation accuracy');
legend(method_names, 'Location', 'southwest', 'NumColumns', 1);

nexttile;
for method = 1:num_methods
    plot(snr_dB, average_selected_pilots(:, method), styles{method}, ...
        'Color', method_colors{method});
    hold on;
end
grid on;
box on;
xlabel('Training SNR (dB)');
ylabel('Average pilot symbols');
title('(b) Training overhead');
ylim([cfg.min_cycles * cfg.num_reference_states - 2, ...
    cfg.max_cycles * cfg.num_reference_states + 2]);
legend(method_names, 'Location', 'east', 'NumColumns', 1);

export_publication_figure(fig, result_dir, ...
    'fig_spatial_stability_guided_training');

fig = publication_figure([80, 80, 940, 350]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
semilogy(snr_dB, max(end_to_end_ser(:,1), 1e-5), '^:', ...
    'Color', colors.gray); hold on;
semilogy(snr_dB, max(end_to_end_ser(:,2), 1e-5), 'o-', ...
    'Color', colors.blue);
semilogy(snr_dB, max(end_to_end_ser(:,3), 1e-5), 's--', ...
    'Color', colors.orange);
semilogy(snr_dB, max(end_to_end_ser(:,4), 1e-5), 'k-.');
grid on; box on; xlabel('Training SNR (dB)'); ylabel('End-to-end ML SER');
title('(a) Detection with estimated CSI');
legend('Fixed short','Stability-guided','Fixed long','Perfect CSI', ...
    'Location','southwest');
nexttile;
for method = [short_method, proposed_method, long_method]
    plot(snr_dB, channel_nmse_dB(:, method), styles{method}, ...
        'Color', method_colors{method}); hold on;
end
grid on; box on; xlabel('Training SNR (dB)'); ylabel('Channel NMSE (dB)');
title('(b) Channel-estimation accuracy');
legend('Fixed short','Stability-guided','Fixed long','Location','southwest');
export_publication_figure(fig, result_dir, 'fig_spatial_end_to_end');

writetable(table(snr_dB(:),end_to_end_ser(:,1),end_to_end_ser(:,2), ...
    end_to_end_ser(:,3),end_to_end_ser(:,4), 'VariableNames', ...
    {'Training_SNR_dB','Fixed_short_SER','Stability_guided_SER', ...
    'Fixed_long_SER','Perfect_CSI_SER'}), ...
    fullfile(result_dir,'spatial_end_to_end_summary.csv'));

save(fullfile(result_dir, ...
    'stability_guided_training_baseline_results.mat'), ...
    'cfg', 'snr_dB', 'method_names', 'channel_nmse_dB', ...
    'channel_nmse_ci95_dB', 'doa_rmse_deg', ...
    'success_probability', 'average_selected_pilots', ...
    'median_selected_pilots', 'selection_histogram', ...
    'average_predicted_channel_nmse_dB', 'end_to_end_ser', 'data_snr_dB');

fprintf('\nAll outputs saved to:\n%s\n', result_dir);

%% Local functions
function pilot_cycles = make_nested_qpsk_pilots(cfg)
    pattern = repmat([1; -1; 1j; -1j], ...
        ceil((cfg.max_cycles * cfg.num_reference_states) / 4), 1);
    pattern = pattern(1:cfg.max_cycles * cfg.num_reference_states);
    first_pilot = ones(cfg.max_cycles * cfg.num_reference_states, 1);
    second_pilot = first_pilot .* pattern;
    pilot_cycles = [first_pilot, second_pilot];
    for num_cycles = cfg.min_cycles:cfg.max_cycles
        num_pilots = num_cycles * cfg.num_reference_states;
        assert(rank(pilot_cycles(1:num_pilots, :)) == cfg.num_users, ...
            'Every candidate pilot prefix must have full column rank.');
    end
end

function S = symbol_codebook(alphabet, L)
    grids = cell(1, L);
    [grids{:}] = ndgrid(alphabet);
    S = zeros(L, numel(grids{1}));
    for user = 1:L
        S(user, :) = grids{user}(:).';
    end
end

function [S, B] = build_spatial_reference_training(pilot_cycles, cfg)
    S = pilot_cycles;
    spatial_phase = (pi / 2) * mod(0:cfg.num_ris - 1, 4);
    B = repmat(cfg.reference_amplitude * exp(1j * spatial_phase), ...
        size(S, 1), 1);
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
    spatial_phase = (pi / 2) * mod(0:cfg.num_ris - 1, 4);
    reference = cfg.reference_amplitude * exp(1j * spatial_phase);
    H_unstructured = zeros(cfg.num_ris, cfg.num_users);
    for element = 1:cfg.num_ris
        weighted_pilots = conj(reference(element)) * pilot_cycles;
        A = 2 * [real(weighted_pilots), -imag(weighted_pilots)];
        y = Z(:, element) - abs(reference(element))^2;
        h_real = (A' * A + 1e-5 * eye(2 * cfg.num_users)) \ (A' * y);
        H_unstructured(element, :) = ...
            (h_real(1:cfg.num_users) + ...
            1j * h_real(cfg.num_users + 1:end)).';
    end
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
    num_parameters = numel(eta);
    degrees_of_freedom = max(numel(z) - num_parameters, 1);
    information.iterations = iteration;
    information.rejections = rejections;
    information.reduced_residual = norm(residual)^2 / ...
        max(degrees_of_freedom * cfg.current_noise_variance, eps);
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
