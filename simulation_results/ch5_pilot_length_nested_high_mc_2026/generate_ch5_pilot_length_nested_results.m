clear;
clc;
close all;

%% Output directory
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
result_dir = fullfile(script_dir, 'results');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

%% Fixed simulation configuration
rng(20260902, 'twister');

cfg.num_users = 2;
cfg.ris_rows = 8;
cfg.ris_cols = 8;
cfg.num_ris = cfg.ris_rows * cfg.ris_cols;
cfg.lambda = 1;
cfg.spacing = cfg.lambda / 2;
cfg.reference_amplitude = 1.5;
cfg.fixed_snr_dB = 12;
cfg.num_monte_carlo = 300;

cfg.theta_true_deg = [15; 35];
cfg.phi_true_deg = [8; -10];
cfg.path_gains = [exp(1j * 0.35); 0.85 * exp(-1j * 0.55)];

cfg.gn_max_iter = 45;
cfg.gn_tolerance = 1e-7;
cfg.gn_initial_damping = 1e-2;
cfg.gn_num_restarts = 2;

pilot_length_vec = [4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 40, 48];
max_pilots = max(pilot_length_vec);

%% Nested pilots and nested reference measurements
% Every operating point uses a prefix of the same maximum-length pilot
% sequence. The reference matrix is likewise designed once and truncated,
% so changing T_p only appends measurements rather than changing their type.
pilot_symbols_full = generate_qpsk_pilots( ...
    max_pilots, cfg.num_users, 20260903);
H_true = generate_channel_matrix(cfg);

reference_phase_full = select_nested_reference_phase( ...
    pilot_symbols_full, H_true, pilot_length_vec, ...
    cfg.reference_amplitude, 128, 20260904);
reference_full = cfg.reference_amplitude * reference_phase_full;

num_lengths = numel(pilot_length_vec);
S_cells = cell(num_lengths, 1);
B_cells = cell(num_lengths, 1);
mu_cells = cell(num_lengths, 1);
noise_variance = zeros(num_lengths, 1);
jacobian_sigma_min = zeros(num_lengths, 1);

for i_length = 1:num_lengths
    num_pilots = pilot_length_vec(i_length);
    S_cells{i_length} = pilot_symbols_full(1:num_pilots, :);
    B_cells{i_length} = reference_full(1:num_pilots, :);
    mu_cells{i_length} = abs( ...
        S_cells{i_length} * H_true.' + B_cells{i_length}).^2;
    noise_variance(i_length) = mean(mu_cells{i_length}(:).^2) / ...
        10^(cfg.fixed_snr_dB / 10);
    [jacobian_sigma_min(i_length), ~] = channel_jacobian_metrics( ...
        S_cells{i_length}, H_true, B_cells{i_length});
end

%% Paired high-Monte-Carlo pilot-length sweep
% The same maximum-length standard-normal noise realization is truncated at
% every T_p. The randomized GN restarts also use the same seed across pilot
% lengths within each trial. Pairing reduces comparison variance without
% changing the marginal distribution at any operating point.
nmse_trials = zeros(cfg.num_monte_carlo, num_lengths);

fprintf('\nNested pilot-length simulation\n');
fprintf('RIS: %d x %d, users: %d, Monte Carlo trials: %d\n', ...
    cfg.ris_rows, cfg.ris_cols, cfg.num_users, cfg.num_monte_carlo);

for mc = 1:cfg.num_monte_carlo
    rng(720000 + mc, 'twister');
    standardized_noise_full = randn(max_pilots, cfg.num_ris);

    for i_length = 1:num_lengths
        num_pilots = pilot_length_vec(i_length);
        z_obs = mu_cells{i_length} + sqrt(noise_variance(i_length)) * ...
            standardized_noise_full(1:num_pilots, :);

        rng(730000 + mc, 'twister');
        H_hat = recover_multiuser_channel_intensity_gn( ...
            z_obs, S_cells{i_length}, B_cells{i_length}, cfg);
        nmse_trials(mc, i_length) = norm(H_hat - H_true, 'fro')^2 / ...
            max(norm(H_true, 'fro')^2, eps);
    end

    if mod(mc, 25) == 0 || mc == cfg.num_monte_carlo
        fprintf('Completed %d/%d trials\n', mc, cfg.num_monte_carlo);
    end
end

nmse_mean_linear = mean(nmse_trials, 1).';
nmse_sem_linear = std(nmse_trials, 0, 1).' / ...
    sqrt(cfg.num_monte_carlo);
nmse_ci_lower_linear = max( ...
    nmse_mean_linear - 1.96 * nmse_sem_linear, eps);
nmse_ci_upper_linear = nmse_mean_linear + 1.96 * nmse_sem_linear;

channel_nmse_pilot_dB = 10 * log10(max(nmse_mean_linear, eps));
nmse_ci_lower_dB = 10 * log10(nmse_ci_lower_linear);
nmse_ci_upper_dB = 10 * log10(nmse_ci_upper_linear);

pilot_table = table( ...
    pilot_length_vec(:), channel_nmse_pilot_dB, ...
    nmse_ci_lower_dB, nmse_ci_upper_dB, ...
    jacobian_sigma_min, noise_variance, ...
    'VariableNames', {'PilotLength', 'Channel_NMSE_dB', ...
    'NMSE_95CI_Lower_dB', 'NMSE_95CI_Upper_dB', ...
    'JacobianSigmaMin', 'NoiseVariance'});
writetable(pilot_table, fullfile(result_dir, ...
    'pilot_length_nested_high_mc.csv'));

fig = publication_figure([100, 100, 540, 390]);
plot(pilot_length_vec, channel_nmse_pilot_dB, 'o-', ...
    'Color', [0, 0.4470, 0.7410], ...
    'MarkerFaceColor', 'w', 'LineWidth', 1.6, 'MarkerSize', 5.5);
grid on;
box on;
xlabel('Pilot length, T_p');
ylabel('Channel NMSE (dB)');
xticks([4, 8, 12, 16, 24, 32, 40, 48]);
xlim([pilot_length_vec(1), pilot_length_vec(end)]);
export_publication_figure(fig, result_dir, ...
    'fig_channel_nmse_pilot_length_nested_high_mc');

save(fullfile(result_dir, ...
    'chapter5_pilot_length_nested_high_mc_results.mat'), ...
    'cfg', 'pilot_length_vec', 'channel_nmse_pilot_dB', ...
    'nmse_ci_lower_dB', 'nmse_ci_upper_dB', 'nmse_trials', ...
    'jacobian_sigma_min', 'noise_variance', ...
    'pilot_symbols_full', 'reference_full', 'H_true');

fprintf('\nPilot-length results:\n');
disp(pilot_table);
fprintf('Number of upward NMSE steps: %d\n', ...
    sum(diff(channel_nmse_pilot_dB) > 0));
fprintf('All outputs saved to:\n%s\n', result_dir);

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

function phase_best = select_nested_reference_phase( ...
        S_full, H, pilot_lengths, reference_amplitude, ...
        num_candidates, seed)
    num_pilots = size(S_full, 1);
    num_ris = size(H, 1);
    phase_best = make_four_phase_reference(num_pilots, num_ris);
    best_score = nested_reference_score( ...
        phase_best, S_full, H, pilot_lengths, reference_amplitude);

    for candidate = 1:num_candidates
        phase_candidate = make_random_reference_phase( ...
            num_pilots, num_ris, seed + candidate);
        score = nested_reference_score( ...
            phase_candidate, S_full, H, pilot_lengths, ...
            reference_amplitude);
        if score > best_score
            best_score = score;
            phase_best = phase_candidate;
        end
    end
end

function score = nested_reference_score( ...
        phase_reference, S_full, H, pilot_lengths, reference_amplitude)
    normalized_sigma_min = zeros(numel(pilot_lengths), 1);
    for i_length = 1:numel(pilot_lengths)
        num_pilots = pilot_lengths(i_length);
        S = S_full(1:num_pilots, :);
        B = reference_amplitude * ...
            phase_reference(1:num_pilots, :);
        sigma_min = channel_jacobian_metrics(S, H, B);
        normalized_sigma_min(i_length) = sigma_min / sqrt(num_pilots);
    end
    score = min(normalized_sigma_min);
end

function phase_reference = make_random_reference_phase( ...
        num_pilots, num_ris, seed)
    old_rng = rng;
    rng(seed, 'twister');
    phase_reference = exp(1j * 2 * pi * rand(num_pilots, num_ris));
    rng(old_rng);
end

function phase_reference = make_four_phase_reference(num_pilots, num_ris)
    time_phase = (pi / 2) * mod((0:num_pilots - 1).', 4);
    space_phase = 2 * pi * (0:num_ris - 1) / num_ris;
    phase_reference = exp(1j * (time_phase + space_phase));
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
    regularization = 1e-3 * trace(B_real.' * B_real) / ...
        max(2 * num_users, 1);

    if norm(B_real, 'fro') > 1e-12
        x_linear = (B_real.' * B_real + ...
            (regularization + 1e-9) * eye(2 * num_users)) \ ...
            (B_real.' * linear_target);
        h_linear = x_linear(1:num_users) + ...
            1j * x_linear(num_users + 1:end);
    else
        h_linear = zeros(num_users, 1);
    end

    signal_scale = sqrt(max(mean(max(z, 0)), 1e-3) / ...
        max(num_users, 1));
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
                    damping * diagonal_scaling + ...
                    1e-10 * eye(2 * num_users);
                delta_candidate = -system_matrix \ gradient;
                h_candidate = h_current + ...
                    delta_candidate(1:num_users) + ...
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
