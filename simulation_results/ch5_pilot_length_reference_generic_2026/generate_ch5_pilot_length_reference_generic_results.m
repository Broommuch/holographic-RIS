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

%% Configuration
% The counting conditions are per RIS element.  We therefore simulate the
% atomic two-user affine phase-retrieval problem y=|X_p h+b|^2.  The full RIS
% model consists of independent copies of this problem with a common pilot
% matrix, so the thresholds are unchanged while substantially more Monte
% Carlo trials and multistart runs can be used here.
cfg.num_users = 2;
cfg.reference_amplitude = 1.5;
cfg.fixed_snr_dB = 12;
cfg.num_monte_carlo = 10000;
cfg.noiseless_restarts = 10;
cfg.noisy_restarts = 3;
cfg.gn_max_iter = 60;
cfg.gn_tolerance = 1e-9;
cfg.gn_initial_damping = 1e-2;
cfg.success_tolerance = 1e-3;

pilot_length_vec = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48];
max_pilots = max(pilot_length_vec);

reference_names = {'No reference', 'Constant', 'Random phase', ...
    'Four-phase'};
num_designs = numel(reference_names);
num_lengths = numel(pilot_length_vec);

%% One nested pilot sequence for every operating point
% The first two rows are deliberately linearly independent.  Every T_p then
% uses a prefix of the same QPSK sequence, so pilot length is the only pilot
% variable in the sweep.
pilot_symbols_full = generate_nested_qpsk_pilots( ...
    max_pilots, cfg.num_users, 20260918);

success_trials = false(cfg.num_monte_carlo, num_lengths, num_designs);
nmse_trials = zeros(cfg.num_monte_carlo, num_lengths, num_designs);
full_rank_trials = false(cfg.num_monte_carlo, num_lengths, num_designs);
jacobian_sigma_min_trials = zeros( ...
    cfg.num_monte_carlo, num_lengths, num_designs);

fprintf('\nPilot-count and reference-diversity simulation\n');
fprintf('Users: %d, trials: %d, SNR: %.1f dB\n', ...
    cfg.num_users, cfg.num_monte_carlo, cfg.fixed_snr_dB);

%% Paired Monte Carlo sweep
for mc = 1:cfg.num_monte_carlo
    rng(810000 + mc, 'twister');
    h_true = (randn(cfg.num_users, 1) + ...
        1j * randn(cfg.num_users, 1)) / sqrt(2);
    h_true = sqrt(cfg.num_users) * h_true / max(norm(h_true), eps);
    standardized_noise_full = randn(max_pilots, 1);
    random_reference_full = cfg.reference_amplitude * ...
        exp(1j * 2 * pi * rand(max_pilots, 1));

    reference_full = cell(num_designs, 1);
    reference_full{1} = zeros(max_pilots, 1);
    reference_full{2} = cfg.reference_amplitude * ones(max_pilots, 1);
    reference_full{3} = random_reference_full;
    reference_full{4} = cfg.reference_amplitude * ...
        exp(1j * (pi / 2) * mod((0:max_pilots - 1).', 4));

    for i_length = 1:num_lengths
        num_pilots = pilot_length_vec(i_length);
        S = pilot_symbols_full(1:num_pilots, :);

        for i_design = 1:num_designs
            b = reference_full{i_design}(1:num_pilots);
            mu = abs(S * h_true + b).^2;

            [jacobian_sigma_min_trials(mc, i_length, i_design), ...
                full_rank_trials(mc, i_length, i_design)] = ...
                channel_jacobian_metrics(S, h_true, b);

            rng(820000 + 10000 * mc + 100 * i_design + i_length, ...
                'twister');
            h_hat_noiseless = recover_affine_intensity_gn( ...
                mu, S, b, cfg, cfg.noiseless_restarts);
            relative_error = norm(h_hat_noiseless - h_true) / ...
                max(norm(h_true), eps);
            success_trials(mc, i_length, i_design) = ...
                relative_error <= cfg.success_tolerance;

            noise_variance = mean(mu.^2) / ...
                10^(cfg.fixed_snr_dB / 10);
            z_noisy = mu + sqrt(noise_variance) * ...
                standardized_noise_full(1:num_pilots);

            rng(830000 + 10000 * mc + 100 * i_design + i_length, ...
                'twister');
            h_hat_noisy = recover_affine_intensity_gn( ...
                z_noisy, S, b, cfg, cfg.noisy_restarts);
            nmse_trials(mc, i_length, i_design) = ...
                norm(h_hat_noisy - h_true)^2 / ...
                max(norm(h_true)^2, eps);
        end
    end

    if mod(mc, 25) == 0 || mc == cfg.num_monte_carlo
        fprintf('Completed %d/%d trials\n', mc, cfg.num_monte_carlo);
    end
end

%% Aggregate and save numerical results
recovery_probability = squeeze(mean(success_trials, 1));
nmse_mean_linear = squeeze(mean(nmse_trials, 1));
channel_nmse_dB = 10 * log10(max(nmse_mean_linear, eps));
full_rank_probability = squeeze(mean(full_rank_trials, 1));
median_jacobian_sigma_min = squeeze( ...
    median(jacobian_sigma_min_trials, 1));

nmse_sem_linear = squeeze(std(nmse_trials, 0, 1)) / ...
    sqrt(cfg.num_monte_carlo);
nmse_ci_lower_dB = 10 * log10(max( ...
    nmse_mean_linear - 1.96 * nmse_sem_linear, eps));
nmse_ci_upper_dB = 10 * log10(max( ...
    nmse_mean_linear + 1.96 * nmse_sem_linear, eps));

length_grid = repmat(pilot_length_vec(:), 1, num_designs);
design_grid = repmat(1:num_designs, num_lengths, 1);
summary_table = table( ...
    length_grid(:), ...
    reshape(string(reference_names(design_grid(:))), [], 1), ...
    recovery_probability(:), channel_nmse_dB(:), ...
    nmse_ci_lower_dB(:), nmse_ci_upper_dB(:), ...
    full_rank_probability(:), median_jacobian_sigma_min(:), ...
    'VariableNames', {'PilotLength', 'ReferenceDesign', ...
    'NoiselessRecoveryProbability', 'ChannelNMSEdB', ...
    'NMSE95CILowerdB', 'NMSE95CIUpperdB', ...
    'JacobianFullRankProbability', 'MedianJacobianSigmaMin'});
writetable(summary_table, fullfile(result_dir, ...
    'pilot_length_reference_generic_summary.csv'));

%% Horizontal two-panel publication figure
colors = [ ...
    0.20, 0.20, 0.20; ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980; ...
    0.4660, 0.6740, 0.1880];
line_styles = {':', '--', '-.', '-'};
markers = {'x', 'o', 'd', 's'};
thresholds = [cfg.num_users, 2 * cfg.num_users, 3 * cfg.num_users];
threshold_labels = {'L', '2L', '3L'};

fig = publication_figure([100, 100, 780, 330]);
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax1 = nexttile(layout, 1);
hold(ax1, 'on');
for i_design = 1:num_designs
    plot(ax1, pilot_length_vec, recovery_probability(:, i_design), ...
        'Color', colors(i_design, :), ...
        'LineStyle', line_styles{i_design}, ...
        'Marker', markers{i_design}, ...
        'LineWidth', 1.45, 'MarkerSize', 4.8, ...
        'MarkerFaceColor', 'w');
end
add_counting_thresholds(ax1, thresholds, threshold_labels);
set(ax1, 'XScale', 'log');
grid(ax1, 'on');
box(ax1, 'on');
xlabel(ax1, 'Pilot length, T_p');
ylabel(ax1, 'Recovery probability');
title(ax1, '(a) Noiseless recovery');
xlim(ax1, [pilot_length_vec(1), pilot_length_vec(end)]);
ylim(ax1, [0, 1.04]);
xticks(ax1, [1, 2, 4, 6, 12, 24, 48]);
xticklabels(ax1, {'1', '2', '4', '6', '12', '24', '48'});
legend(ax1, reference_names, 'Location', 'southeast', ...
    'FontSize', 7.5);

ax2 = nexttile(layout, 2);
hold(ax2, 'on');
for i_design = 1:num_designs
    plot(ax2, pilot_length_vec, channel_nmse_dB(:, i_design), ...
        'Color', colors(i_design, :), ...
        'LineStyle', line_styles{i_design}, ...
        'Marker', markers{i_design}, ...
        'LineWidth', 1.45, 'MarkerSize', 4.8, ...
        'MarkerFaceColor', 'w');
end
add_counting_thresholds(ax2, thresholds, threshold_labels);
set(ax2, 'XScale', 'log');
grid(ax2, 'on');
box(ax2, 'on');
xlabel(ax2, 'Pilot length, T_p');
ylabel(ax2, 'Channel NMSE (dB)');
title(ax2, '(b) SNR = 12 dB');
xlim(ax2, [pilot_length_vec(1), pilot_length_vec(end)]);
xticks(ax2, [1, 2, 4, 6, 12, 24, 48]);
xticklabels(ax2, {'1', '2', '4', '6', '12', '24', '48'});

export_publication_figure(fig, result_dir, ...
    'fig_pilot_length_reference_generic');

save(fullfile(result_dir, ...
    'chapter5_pilot_length_reference_generic_results.mat'), ...
    'cfg', 'pilot_length_vec', 'reference_names', ...
    'pilot_symbols_full', 'recovery_probability', ...
    'channel_nmse_dB', 'nmse_ci_lower_dB', 'nmse_ci_upper_dB', ...
    'full_rank_probability', 'median_jacobian_sigma_min', ...
    'success_trials', 'nmse_trials', 'jacobian_sigma_min_trials');

fprintf('\nNoiseless recovery probability:\n');
disp(array2table(recovery_probability, ...
    'VariableNames', matlab.lang.makeValidName(reference_names), ...
    'RowNames', compose('Tp_%d', pilot_length_vec)));
fprintf('\nChannel NMSE at %.1f dB:\n', cfg.fixed_snr_dB);
disp(array2table(channel_nmse_dB, ...
    'VariableNames', matlab.lang.makeValidName(reference_names), ...
    'RowNames', compose('Tp_%d', pilot_length_vec)));
fprintf('All outputs saved to:\n%s\n', result_dir);

%% Local functions
function S = generate_nested_qpsk_pilots(num_pilots, num_users, seed)
    old_rng = rng;
    rng(seed, 'twister');
    constellation = [1 + 1j; -1 + 1j; -1 - 1j; 1 - 1j] / sqrt(2);
    indices = randi([1, 4], num_pilots, num_users);
    S = constellation(indices);

    if num_users == 2 && num_pilots >= 2
        S(1, :) = [constellation(1), constellation(1)];
        S(2, :) = [constellation(1), constellation(2)];
    end
    rng(old_rng);
end

function [sigma_min, is_full_rank] = channel_jacobian_metrics(S, h, b)
    num_real_parameters = 2 * size(S, 2);
    q = S * h + b;
    weighted_S = bsxfun(@times, conj(q), S);
    J = 2 * [real(weighted_S), -imag(weighted_S)];
    singular_values = svd(J, 'econ');
    rank_tolerance = max(size(J)) * eps(max(singular_values));
    numerical_rank = sum(singular_values > rank_tolerance);
    is_full_rank = numerical_rank == num_real_parameters;

    % An underdetermined Jacobian has an unreported zero singular value when
    % svd(...,'econ') is used.  Record that zero explicitly.
    if is_full_rank
        sigma_min = singular_values(end);
    else
        sigma_min = 0;
    end
end

function h_best = recover_affine_intensity_gn( ...
        z, S, b, cfg, num_restarts)
    num_users = size(S, 2);
    z = z(:);
    b = b(:);

    weighted_S = bsxfun(@times, conj(b), S);
    A_linear = 2 * [real(weighted_S), -imag(weighted_S)];
    linear_target = z - abs(b).^2;
    linear_regularization = 1e-4 * trace(A_linear.' * A_linear) / ...
        max(2 * num_users, 1);

    if norm(A_linear, 'fro') > 1e-12
        x_linear = (A_linear.' * A_linear + ...
            (linear_regularization + 1e-10) * eye(2 * num_users)) \ ...
            (A_linear.' * linear_target);
        h_linear = x_linear(1:num_users) + ...
            1j * x_linear(num_users + 1:end);
    else
        h_linear = zeros(num_users, 1);
    end

    signal_scale = sqrt(max(mean(max(z, 0)), 1e-6) / ...
        max(num_users, 1));
    h_best = h_linear;
    cost_best = inf;

    for restart = 1:num_restarts
        if restart == 1 && norm(h_linear) > 1e-12
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

            for damping_trial = 1:12
                system_matrix = normal_matrix + ...
                    damping * diagonal_scaling + ...
                    1e-11 * eye(2 * num_users);
                delta_candidate = -system_matrix \ gradient;
                h_candidate = h_current + ...
                    delta_candidate(1:num_users) + ...
                    1j * delta_candidate(num_users + 1:end);
                cost_candidate = intensity_cost(h_candidate, S, b, z);

                if cost_candidate < cost_current
                    h_current = h_candidate;
                    cost_current = cost_candidate;
                    delta = delta_candidate;
                    damping = max(damping / 3, 1e-12);
                    accepted = true;
                    break;
                end
                damping = min(damping * 10, 1e12);
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

function add_counting_thresholds(ax, thresholds, labels)
    styles = {':', '--', '-.'};
    for i = 1:numel(thresholds)
        xline(ax, thresholds(i), styles{i}, labels{i}, ...
            'Color', [0.45, 0.45, 0.45], ...
            'LineWidth', 0.8, 'LabelVerticalAlignment', 'top', ...
            'LabelHorizontalAlignment', 'left', ...
            'LabelOrientation', 'horizontal', ...
            'HandleVisibility', 'off');
    end
end

function fig = publication_figure(position)
    fig = figure('Color', 'w', 'Position', position, 'Visible', 'off');
    set(groot, 'defaultAxesFontName', 'Times New Roman');
    set(groot, 'defaultAxesFontSize', 9);
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
