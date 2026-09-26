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

%% Configuration inherited from the original reference-design experiment
rng(20260903, 'twister');

cfg.num_users = 2;
cfg.ris_rows = 8;
cfg.ris_cols = 8;
cfg.num_ris = cfg.ris_rows * cfg.ris_cols;
cfg.lambda = 1;
cfg.spacing = cfg.lambda / 2;
cfg.num_pilots = 32;
cfg.reference_amplitude = 1.5;
cfg.num_monte_carlo = 300;
cfg.fixed_snr_dB = 12;

cfg.theta_true_deg = [15; 35];
cfg.phi_true_deg = [8; -10];
cfg.path_gains = [exp(1j * 0.35); 0.85 * exp(-1j * 0.55)];

cfg.gn_max_iter = 45;
cfg.gn_tolerance = 1e-7;
cfg.gn_initial_damping = 1e-2;
cfg.gn_num_restarts = 2;
cfg.reference_codebook_size = 32;

reference_names = {'No reference', 'Constant', 'Random', ...
    'Uniform', 'Optimized'};
reference_variable_names = {'NoReference', 'Constant', 'Random', ...
    'Uniform', 'Optimized'};

% Preserve the original script's random-number ordering: it first generates
% the 48-row master pilot matrix used by all sweeps and then takes 32 rows.
pilot_symbols_full = generate_qpsk_pilots( ...
    48, cfg.num_users, 20260829);
pilot_symbols = pilot_symbols_full(1:cfg.num_pilots, :);
H_true = generate_channel_matrix(cfg);

reference_phase_optimized = select_optimized_reference_phase( ...
    pilot_symbols, H_true, cfg.reference_codebook_size, 20260830);
reference_baseline = cfg.reference_amplitude * reference_phase_optimized;

num_designs = numel(reference_names);
reference_cells = cell(num_designs, 1);
reference_cells{1} = zeros(cfg.num_pilots, cfg.num_ris);
reference_cells{2} = cfg.reference_amplitude * ...
    ones(cfg.num_pilots, cfg.num_ris);
reference_cells{3} = cfg.reference_amplitude * ...
    make_random_reference_phase( ...
    cfg.num_pilots, cfg.num_ris, 20260831);
reference_cells{4} = cfg.reference_amplitude * ...
    make_uniform_reference_phase(cfg.num_pilots, cfg.num_ris);
reference_cells{5} = reference_baseline;

%% Precompute the noiseless measurements and Jacobian metrics
mu_cells = cell(num_designs, 1);
noise_variance = zeros(num_designs, 1);
jacobian_sigma_min = zeros(num_designs, 1);
jacobian_condition_number = zeros(num_designs, 1);

for i_design = 1:num_designs
    B = reference_cells{i_design};
    mu_cells{i_design} = abs(pilot_symbols * H_true.' + B).^2;
    noise_variance(i_design) = mean(mu_cells{i_design}(:).^2) / ...
        10^(cfg.fixed_snr_dB / 10);
    [jacobian_sigma_min(i_design), ...
        jacobian_condition_number(i_design)] = ...
        channel_jacobian_metrics(pilot_symbols, H_true, B);
end

%% Paired Monte Carlo comparison
% A common standard-normal noise realization and common randomized GN restart
% are used across the five designs in each trial. This preserves the marginal
% distribution for every design while reducing comparison variance.
nmse_trials = zeros(cfg.num_monte_carlo, num_designs);

fprintf('\nReference-design point-estimate simulation\n');
fprintf('RIS: %d x %d, users: %d, pilots: %d, Monte Carlo: %d\n', ...
    cfg.ris_rows, cfg.ris_cols, cfg.num_users, cfg.num_pilots, ...
    cfg.num_monte_carlo);

for mc = 1:cfg.num_monte_carlo
    rng(810000 + mc, 'twister');
    standardized_noise = randn(cfg.num_pilots, cfg.num_ris);

    for i_design = 1:num_designs
        z_obs = mu_cells{i_design} + sqrt(noise_variance(i_design)) * ...
            standardized_noise;

        rng(820000 + mc, 'twister');
        H_hat = recover_multiuser_channel_intensity_gn( ...
            z_obs, pilot_symbols, reference_cells{i_design}, cfg);
        nmse_trials(mc, i_design) = norm(H_hat - H_true, 'fro')^2 / ...
            max(norm(H_true, 'fro')^2, eps);
    end

    if mod(mc, 25) == 0 || mc == cfg.num_monte_carlo
        fprintf('Completed %d/%d trials\n', mc, cfg.num_monte_carlo);
    end
end

%% Point estimates and 95% Monte Carlo confidence intervals
nmse_mean_linear = mean(nmse_trials, 1).';
nmse_sem_linear = std(nmse_trials, 0, 1).' / ...
    sqrt(cfg.num_monte_carlo);
nmse_ci_lower_linear = max( ...
    nmse_mean_linear - 1.96 * nmse_sem_linear, eps);
nmse_ci_upper_linear = nmse_mean_linear + 1.96 * nmse_sem_linear;

channel_nmse_dB = 10 * log10(max(nmse_mean_linear, eps));
nmse_ci_lower_dB = 10 * log10(nmse_ci_lower_linear);
nmse_ci_upper_dB = 10 * log10(nmse_ci_upper_linear);

summary_table = table( ...
    string(reference_names(:)), channel_nmse_dB, ...
    nmse_ci_lower_dB, nmse_ci_upper_dB, ...
    jacobian_sigma_min, jacobian_condition_number, noise_variance, ...
    'VariableNames', {'ReferenceDesign', 'Channel_NMSE_dB', ...
    'NMSE_95CI_Lower_dB', 'NMSE_95CI_Upper_dB', ...
    'JacobianSigmaMin', 'JacobianConditionNumber', 'NoiseVariance'});
writetable(summary_table, fullfile(result_dir, ...
    'reference_design_point_ci_summary.csv'));

trial_table = array2table(nmse_trials, ...
    'VariableNames', reference_variable_names);
writetable(trial_table, fullfile(result_dir, ...
    'reference_design_nmse_trials_linear.csv'));

%% Two-panel point-estimate figure
fig = publication_figure([100, 100, 570, 600]);
layout = tiledlayout(2, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
x = 1:num_designs;

ax_top = nexttile(layout);
lower_error = channel_nmse_dB - nmse_ci_lower_dB;
upper_error = nmse_ci_upper_dB - channel_nmse_dB;
errorbar(ax_top, x, channel_nmse_dB, lower_error, upper_error, ...
    'o', 'LineStyle', 'none', 'Color', [0, 0.4470, 0.7410], ...
    'MarkerFaceColor', 'w', 'LineWidth', 1.5, ...
    'MarkerSize', 6, 'CapSize', 6);
yline(ax_top, 0, ':', 'Color', [0.35, 0.35, 0.35], ...
    'LineWidth', 0.9);
grid(ax_top, 'on');
box(ax_top, 'on');
ylabel(ax_top, 'Channel NMSE (dB)');
set(ax_top, 'XTick', x, 'XTickLabel', [], ...
    'XLim', [0.6, num_designs + 0.4]);
text(ax_top, 0.02, 0.93, '(a)', 'Units', 'normalized', ...
    'FontWeight', 'bold');

ax_bottom = nexttile(layout);
sigma_min_plot = jacobian_sigma_min;
sigma_min_plot(sigma_min_plot < 1e-12) = 0;
plot(ax_bottom, x, sigma_min_plot, ...
    's', 'LineStyle', 'none', 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor', 'w', 'LineWidth', 1.5, 'MarkerSize', 6);
grid(ax_bottom, 'on');
box(ax_bottom, 'on');
ylabel(ax_bottom, 'Minimum Jacobian singular value');
set(ax_bottom, 'XTick', x, 'XTickLabel', reference_names, ...
    'XTickLabelRotation', 18, ...
    'XLim', [0.6, num_designs + 0.4], 'YLim', [-0.5, 10.8]);
text(ax_bottom, 0.02, 0.93, '(b)', 'Units', 'normalized', ...
    'FontWeight', 'bold');

export_publication_figure(fig, result_dir, ...
    'fig_reference_design_point_ci');

save(fullfile(result_dir, ...
    'chapter5_reference_design_point_ci_results.mat'), ...
    'cfg', 'reference_names', 'reference_cells', ...
    'channel_nmse_dB', 'nmse_ci_lower_dB', ...
    'nmse_ci_upper_dB', 'nmse_trials', ...
    'jacobian_sigma_min', 'jacobian_condition_number', ...
    'noise_variance', 'pilot_symbols', 'H_true');

fprintf('\nReference-design results:\n');
disp(summary_table);
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

function phase_reference = make_random_reference_phase( ...
        num_pilots, num_ris, seed)
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

function phase_best = select_optimized_reference_phase( ...
        S, H, num_candidates, seed)
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
