%% Perfect-CSI detection with projected-WF baseline
% Extends the spatial-reference detector experiment with a projected
% Wirtinger-flow baseline under the same affine intensity model,
% initialization, iteration budget, and exact-distance stopping certificate.

clc;
clear;
close all;

script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

rng(20261011, 'twister');

Ny = 8;
Nz = 8;
N = Ny * Nz;
L = 2;
theta = [15; 35];
phi = [8; -10];
alpha = [exp(1j * 0.35); 0.85 * exp(-1j * 0.55)];
H = channel_matrix(Ny, Nz, theta, phi, alpha);
qpsk = exp(1j * (pi / 4 + (0:3) * pi / 2)).';
S = symbol_codebook(qpsk, L);
num_candidates = size(S, 2);

phase = (pi / 2) * mod((0:N - 1).', 4);
b = 1.5 * exp(1j * phase);
Mu = abs(H * S + b).^2;
dmat = pairwise_distances(Mu);
dmin = min(dmat(dmat > 0));

snr_dB = -15:2:15;
num_trials = 30000;
if strcmp(getenv('CODEX_FAST_TEST'), '1')
    num_trials = 2000;
end
max_iterations = 25;

ser_ml = zeros(size(snr_dB));
ser_lin = zeros(size(snr_dB));
ser_pwf = zeros(size(snr_dB));
ser_gs = zeros(size(snr_dB));
union_bound = zeros(size(snr_dB));
iterations_pwf = zeros(size(snr_dB));
iterations_gs = zeros(size(snr_dB));
fixed_budget = max_iterations * ones(size(snr_dB));
spectral_norm_squared = norm(H, 2)^2;

fprintf('\nProjected-WF baseline under the spatial affine model\n');
fprintf('Trials per SNR: %d, candidates: %d\n', num_trials, num_candidates);

for snr_index = 1:numel(snr_dB)
    sigma2 = mean(Mu(:).^2) / 10^(snr_dB(snr_index) / 10);
    sigma = sqrt(sigma2);
    union_bound(snr_index) = symbol_union_bound(dmat, S, sigma);

    transmitted_index = randi(num_candidates, 1, num_trials);
    transmitted_symbols = S(:, transmitted_index);
    z = Mu(:, transmitted_index) + sigma * randn(N, num_trials);

    ml_distance = zeros(num_candidates, num_trials);
    for candidate = 1:num_candidates
        ml_distance(candidate, :) = sum((z - Mu(:, candidate)).^2, 1);
    end
    [~, detected_index] = min(ml_distance, [], 1);
    errors_ml = sum(S(:, detected_index) ~= transmitted_symbols, 'all');

    [s_linear, s_continuous] = linearized_estimate_batch( ...
        z, H, b, qpsk);
    errors_lin = sum(s_linear ~= transmitted_symbols, 'all');

    [s_pwf, used_pwf] = projected_wf_detect_batch( ...
        z, H, b, qpsk, s_continuous, max_iterations, dmin, ...
        spectral_norm_squared);
    errors_pwf = sum(s_pwf ~= transmitted_symbols, 'all');

    [s_gs, used_gs] = gs_detect_batch( ...
        z, H, b, qpsk, s_linear, max_iterations, dmin);
    errors_gs = sum(s_gs ~= transmitted_symbols, 'all');

    ser_ml(snr_index) = errors_ml / (num_trials * L);
    ser_lin(snr_index) = errors_lin / (num_trials * L);
    ser_pwf(snr_index) = errors_pwf / (num_trials * L);
    ser_gs(snr_index) = errors_gs / (num_trials * L);
    iterations_pwf(snr_index) = mean(used_pwf);
    iterations_gs(snr_index) = mean(used_gs);

    fprintf(['SNR %5.1f dB: ML %.4g, linearized %.4g, PWF %.4g, ' ...
        'GS %.4g, iterations %.2f/%.2f\n'], snr_dB(snr_index), ...
        ser_ml(snr_index), ser_lin(snr_index), ser_pwf(snr_index), ...
        ser_gs(snr_index), iterations_pwf(snr_index), ...
        iterations_gs(snr_index));
end

% Zero empirical SER cannot be represented on a logarithmic axis. Omit those
% points instead of replacing them by an artificial plotting floor.
ser_ml_plot = ser_ml;
ser_lin_plot = ser_lin;
ser_pwf_plot = ser_pwf;
ser_gs_plot = ser_gs;
ser_ml_plot(ser_ml_plot == 0) = NaN;
ser_lin_plot(ser_lin_plot == 0) = NaN;
ser_pwf_plot(ser_pwf_plot == 0) = NaN;
ser_gs_plot(ser_gs_plot == 0) = NaN;
fig = publication_figure([100, 100, 920, 350]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
semilogy(snr_dB, ser_ml_plot, 'o-', ...
    'Color', [0, 0.447, 0.741]);
hold on;
semilogy(snr_dB, ser_lin_plot, '^-', ...
    'Color', [0.466, 0.674, 0.188]);
semilogy(snr_dB, ser_pwf_plot, 'd-.', ...
    'Color', [0.494, 0.184, 0.556]);
semilogy(snr_dB, ser_gs_plot, 's--', ...
    'Color', [0.85, 0.325, 0.098]);
semilogy(snr_dB, union_bound, 'k-.');
grid on;
box on;
xlabel('SNR (dB)');
ylabel('SER');
ylim([1e-5, 1]);
title('(a) Detection accuracy');
legend('ML', 'Linearized', 'Projected WF', 'Reference-assisted GS', ...
    'Union bound', 'Location', 'southwest', 'FontSize', 8.5);

nexttile;
plot(snr_dB, fixed_budget, 'k--');
hold on;
plot(snr_dB, iterations_pwf, 'd-.', ...
    'Color', [0.494, 0.184, 0.556]);
plot(snr_dB, iterations_gs, 'o-', ...
    'Color', [0, 0.447, 0.741]);
grid on;
box on;
xlabel('SNR (dB)');
ylabel('Average nonlinear iterations');
title('(b) Certified stopping');
legend('Fixed budget', 'Projected WF', 'Reference-assisted GS', ...
    'Location', 'southwest', 'FontSize', 8.5);

export_figure(fig, result_dir, 'fig_projected_wf_detection');

summary = table(snr_dB(:), ser_ml(:), ser_lin(:), ser_pwf(:), ser_gs(:), ...
    union_bound(:), iterations_pwf(:), iterations_gs(:), ...
    'VariableNames', {'SNR_dB', 'ML_SER', 'Linearized_SER', ...
    'Projected_WF_SER', 'GS_SER', 'Union_bound', ...
    'Projected_WF_iterations', 'GS_iterations'});
writetable(summary, fullfile(result_dir, ...
    'projected_wf_detection_summary.csv'));

save(fullfile(result_dir, 'projected_wf_detection_results.mat'), ...
    'snr_dB', 'ser_ml', 'ser_lin', 'ser_pwf', 'ser_gs', ...
    'union_bound', 'iterations_pwf', 'iterations_gs', ...
    'H', 'b', 'S', 'dmin', 'num_trials', 'max_iterations');

function H = channel_matrix(Ny, Nz, theta_deg, phi_deg, alpha)
    [yy, zz] = meshgrid(0:Nz - 1, 0:Ny - 1);
    y = yy(:) - (Nz - 1) / 2;
    z = zz(:) - (Ny - 1) / 2;
    H = zeros(Ny * Nz, numel(theta_deg));
    for user = 1:numel(theta_deg)
        theta = theta_deg(user) * pi / 180;
        phi = phi_deg(user) * pi / 180;
        steering = exp(1j * pi * sin(theta) .* ...
            (y * cos(phi) + z * sin(phi)));
        H(:, user) = alpha(user) * steering;
    end
end

function S = symbol_codebook(alphabet, num_users)
    grids = cell(1, num_users);
    [grids{:}] = ndgrid(alphabet);
    S = zeros(num_users, numel(grids{1}));
    for user = 1:num_users
        S(user, :) = grids{user}(:).';
    end
end

function distances = pairwise_distances(codewords)
    num_candidates = size(codewords, 2);
    distances = inf(num_candidates);
    for first = 1:num_candidates
        for second = first + 1:num_candidates
            value = norm(codewords(:, first) - codewords(:, second));
            distances(first, second) = value;
            distances(second, first) = value;
        end
    end
end

function bound = symbol_union_bound(distances, symbols, sigma)
    [num_users, num_candidates] = size(symbols);
    bound = 0;
    for first = 1:num_candidates
        for second = 1:num_candidates
            if first ~= second
                symbol_fraction = ...
                    sum(symbols(:, first) ~= symbols(:, second)) / num_users;
                bound = bound + symbol_fraction * 0.5 * ...
                    erfc(distances(first, second) / (2 * sigma * sqrt(2)));
            end
        end
    end
    bound = min(bound / num_candidates, 1);
end

function [symbols, continuous] = linearized_estimate_batch(z, G, b, alphabet)
    A = 2 * [real(conj(b) .* G), -imag(conj(b) .* G)];
    operator = (A' * A + 1e-6 * eye(size(A, 2))) \ A';
    estimate = operator * (z - abs(b).^2);
    num_users = size(G, 2);
    continuous = estimate(1:num_users, :) + ...
        1j * estimate(num_users + 1:end, :);
    symbols = constellation_projection(continuous, alphabet);
end

function [symbols, used] = projected_wf_detect_batch( ...
        z, G, b, alphabet, initial, max_iterations, dmin, ...
        spectral_norm_squared)
    continuous = convex_hull_projection(initial, alphabet);
    symbols = constellation_projection(continuous, alphabet);
    num_trials = size(z, 2);
    used = zeros(1, num_trials);
    certificate = sqrt(sum((z - abs(G * symbols + b).^2).^2, 1)) < dmin / 2;
    active = ~certificate;

    for iteration = 1:max_iterations
        active_index = find(active);
        if isempty(active_index)
            break;
        end
        current = continuous(:, active_index);
        z_active = z(:, active_index);
        field = G * current + b;
        residual = abs(field).^2 - z_active;
        cost = 0.5 * sum(residual.^2, 1);
        gradient = 2 * G' * (residual .* field);
        scale = spectral_norm_squared * ...
            (2 * max(abs(field).^2, [], 1) + max(abs(residual), [], 1) + 1e-9);
        step = 0.5 ./ max(scale, 1e-9);

        trial = current;
        accepted = false(1, numel(active_index));
        unresolved = true(1, numel(active_index));
        for line_search = 1:12
            local_index = find(unresolved);
            if isempty(local_index)
                break;
            end
            candidate = convex_hull_projection( ...
                current(:, local_index) - ...
                gradient(:, local_index) .* step(local_index), alphabet);
            candidate_residual = abs(G * candidate + b).^2 - ...
                z_active(:, local_index);
            local_accept = 0.5 * sum(candidate_residual.^2, 1) <= ...
                cost(local_index);
            accepted_index = local_index(local_accept);
            trial(:, accepted_index) = candidate(:, local_accept);
            accepted(accepted_index) = true;
            unresolved(accepted_index) = false;
            rejected_index = local_index(~local_accept);
            step(rejected_index) = step(rejected_index) / 2;
        end

        accepted_global = active_index(accepted);
        stalled_global = active_index(~accepted);
        if ~isempty(stalled_global)
            active(stalled_global) = false;
        end
        if isempty(accepted_global)
            continue;
        end

        relative_change = sqrt(sum(abs(trial(:, accepted) - ...
            current(:, accepted)).^2, 1)) ./ ...
            max(sqrt(sum(abs(current(:, accepted)).^2, 1)), 1e-9);
        continuous(:, accepted_global) = trial(:, accepted);
        symbols(:, accepted_global) = constellation_projection( ...
            continuous(:, accepted_global), alphabet);
        used(accepted_global) = iteration;
        certificate = sqrt(sum((z(:, accepted_global) - ...
            abs(G * symbols(:, accepted_global) + b).^2).^2, 1)) < dmin / 2;
        converged = relative_change < 1e-6;
        active(accepted_global(certificate | converged)) = false;
    end
end

function [symbols, used] = gs_detect_batch( ...
        z, G, b, alphabet, symbols, max_iterations, dmin)
    num_trials = size(z, 2);
    used = zeros(1, num_trials);
    certificate = sqrt(sum((z - abs(G * symbols + b).^2).^2, 1)) < dmin / 2;
    active = ~certificate;
    continuous = symbols;
    back_projection = (G' * G + 1e-6 * eye(size(G, 2))) \ G';

    for iteration = 1:max_iterations
        active_index = find(active);
        if isempty(active_index)
            break;
        end
        field = G * continuous(:, active_index) + b;
        projected_field = sqrt(max(z(:, active_index), 0)) .* ...
            exp(1j * angle(field));
        continuous(:, active_index) = back_projection * ...
            (projected_field - b);
        symbols(:, active_index) = constellation_projection( ...
            continuous(:, active_index), alphabet);
        used(active_index) = iteration;
        certificate = sqrt(sum((z(:, active_index) - ...
            abs(G * symbols(:, active_index) + b).^2).^2, 1)) < dmin / 2;
        active(active_index(certificate)) = false;
    end
end

function projected = convex_hull_projection(values, alphabet)
    real_limit = max(abs(real(alphabet)));
    imag_limit = max(abs(imag(alphabet)));
    projected = min(max(real(values), -real_limit), real_limit) + ...
        1j * min(max(imag(values), -imag_limit), imag_limit);
end

function symbols = constellation_projection(values, alphabet)
    symbols = zeros(size(values));
    for index = 1:numel(values)
        [~, nearest] = min(abs(values(index) - alphabet));
        symbols(index) = alphabet(nearest);
    end
end

function fig = publication_figure(position)
    fig = figure('Color', 'w', 'Position', position, 'Visible', 'off');
    set(groot, 'defaultAxesFontName', 'Times New Roman', ...
        'defaultTextFontName', 'Times New Roman', ...
        'defaultAxesFontSize', 11, 'defaultLineLineWidth', 1.5, ...
        'defaultLineMarkerSize', 6);
end

function export_figure(fig, result_dir, base_name)
    exportgraphics(fig, fullfile(result_dir, [base_name, '.eps']), ...
        'ContentType', 'vector');
    exportgraphics(fig, fullfile(result_dir, [base_name, '.png']), ...
        'Resolution', 300);
    savefig(fig, fullfile(result_dir, [base_name, '.fig']));
    close(fig);
end
