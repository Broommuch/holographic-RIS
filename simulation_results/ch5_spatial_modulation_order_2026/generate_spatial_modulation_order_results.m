%% Modulation-order comparison with the fixed simultaneous spatial reference
clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');
if ~exist(result_dir, 'dir'), mkdir(result_dir); end
rng(20260928, 'twister');

%% Configuration consistent with Chapter V
cfg.num_users = 2;
cfg.ris_rows = 8;
cfg.ris_cols = 8;
cfg.num_ris = cfg.ris_rows * cfg.ris_cols;
cfg.lambda = 1;
cfg.spacing = cfg.lambda / 2;
cfg.reference_amplitude = 1.5;
cfg.theta_true_deg = [15; 35];
cfg.phi_true_deg = [8; -10];
cfg.path_gains = [exp(1j * 0.35); 0.85 * exp(-1j * 0.55)];
cfg.num_trials = 10000;
cfg.trial_chunk_size = 100;
cfg.codebook_chunk_size = 4096;
cfg.basis_sample_size = 1024;
cfg.basis_relative_tolerance = 1e-10;
cfg.projection_residual_tolerance = 1e-9;

snr_dB_vec = -12:3:24;
modulation_orders = [4, 16, 64, 256];
modulation_names = {'QPSK', '16-QAM', '64-QAM', '256-QAM'};

%% Channel and fixed spatial four-phase reference
H_true = generate_channel_matrix(cfg);
spatial_phase = (pi / 2) * mod((0:cfg.num_ris - 1).', 4);
b_spatial = cfg.reference_amplitude * exp(1j * spatial_phase);

num_modulations = numel(modulation_orders);
num_snr = numel(snr_dB_vec);
num_symbol_decisions = cfg.num_trials * cfg.num_users;

ser = zeros(num_snr, num_modulations);
ser_plot = zeros(num_snr, num_modulations);
symbol_error_count = zeros(num_snr, num_modulations);
ser_ci_lower = zeros(num_snr, num_modulations);
ser_ci_upper = zeros(num_snr, num_modulations);
basis_rank = zeros(num_modulations, 1);
projection_residual = zeros(num_modulations, 1);
noiseless_intensity_power = zeros(num_modulations, 1);
num_joint_candidates = zeros(num_modulations, 1);
projected_codebooks = cell(num_modulations, 1);

fprintf('\nSpatial-reference modulation-order simulation\n');
fprintf('Surface: %d x %d, users: %d, trials per SNR: %d\n', ...
    cfg.ris_rows, cfg.ris_cols, cfg.num_users, cfg.num_trials);

%% Exact ML in the lossless intensity-difference subspace
for i_modulation = 1:num_modulations
    order = modulation_orders(i_modulation);
    constellation = generate_square_qam(order);
    [combination_indices, candidates] = enumerate_symbol_vectors( ...
        constellation, cfg.num_users);
    num_joint_candidates(i_modulation) = size(candidates, 1);

    [projected_mu, noiseless_power, rank_current, residual_current] = ...
        build_reduced_intensity_codebook(H_true, b_spatial, candidates, ...
        cfg, 810000 + i_modulation);
    projected_codebooks{i_modulation} = projected_mu;
    noiseless_intensity_power(i_modulation) = noiseless_power;
    basis_rank(i_modulation) = rank_current;
    projection_residual(i_modulation) = residual_current;

    fprintf('%s: %d joint candidates, rank %d, residual %.3e\n', ...
        modulation_names{i_modulation}, num_joint_candidates(i_modulation), ...
        rank_current, residual_current);

    for i_snr = 1:num_snr
        noise_variance = noiseless_power / 10^(snr_dB_vec(i_snr) / 10);
        symbol_error_count(i_snr, i_modulation) = ...
            simulate_exact_projected_ml_errors(projected_mu, ...
            combination_indices, noise_variance, cfg.num_trials, ...
            cfg.trial_chunk_size, 820000 + i_modulation);
        ser(i_snr, i_modulation) = ...
            symbol_error_count(i_snr, i_modulation) / num_symbol_decisions;
        [ser_ci_lower(i_snr, i_modulation), ...
            ser_ci_upper(i_snr, i_modulation)] = wilson_interval( ...
            symbol_error_count(i_snr, i_modulation), ...
            num_symbol_decisions, 1.96);
    end

    ser_plot(:, i_modulation) = ser(:, i_modulation);
    ser_plot(ser_plot(:, i_modulation) == 0, i_modulation) = NaN;
end

%% Save numerical results
[snr_grid, modulation_grid] = ndgrid(snr_dB_vec(:), 1:num_modulations);
modulation_name_column = reshape( ...
    string(modulation_names(modulation_grid(:))), [], 1);
modulation_order_column = reshape( ...
    modulation_orders(modulation_grid(:)), [], 1);
zero_error_indicator = symbol_error_count == 0;

summary_table = table(snr_grid(:), modulation_name_column, ...
    modulation_order_column, num_joint_candidates(modulation_grid(:)), ...
    symbol_error_count(:), ser(:), zero_error_indicator(:), ...
    ser_ci_lower(:), ser_ci_upper(:), basis_rank(modulation_grid(:)), ...
    projection_residual(modulation_grid(:)), ...
    'VariableNames', {'SNR_dB', 'Modulation', 'Order', ...
    'NumberOfJointCandidates', 'SymbolErrorCount', 'EmpiricalSER', ...
    'ZeroErrorPoint', 'SER_95CI_Lower', 'SER_95CI_Upper', ...
    'ReducedSubspaceRank', 'ProjectionResidual'});
writetable(summary_table, fullfile(result_dir, ...
    'spatial_modulation_order_ser_summary.csv'));

save(fullfile(result_dir, 'spatial_modulation_order_results.mat'), ...
    'cfg', 'snr_dB_vec', 'modulation_orders', 'modulation_names', ...
    'ser', 'ser_plot', 'symbol_error_count', 'ser_ci_lower', ...
    'ser_ci_upper', 'zero_error_indicator', 'basis_rank', ...
    'projection_residual', 'noiseless_intensity_power', ...
    'num_joint_candidates', 'projected_codebooks', 'H_true', 'b_spatial');

%% Publication figure: zero-error points are omitted
colors = [0.0000, 0.4470, 0.7410; ...
          0.8500, 0.3250, 0.0980; ...
          0.4660, 0.6740, 0.1880; ...
          0.4940, 0.1840, 0.5560];
line_styles = {'-', '--', '-.', ':'};
markers = {'o', 's', 'd', '^'};

fig = publication_figure([100, 100, 520, 400]);
hold on;
for i_modulation = 1:num_modulations
    semilogy(snr_dB_vec, ser_plot(:, i_modulation), ...
        'Color', colors(i_modulation, :), ...
        'LineStyle', line_styles{i_modulation}, ...
        'Marker', markers{i_modulation}, 'MarkerFaceColor', 'w');
end
set(gca, 'YScale', 'log');
grid on; box on;
xlabel('SNR (dB)');
ylabel('ML SER');
xlim([snr_dB_vec(1), snr_dB_vec(end)]);
ylim([4e-5, 1]);
xticks(-12:6:24);
legend(modulation_names, 'Location', 'southwest', 'FontSize', 9);
export_publication_figure(fig, result_dir, 'fig_spatial_modulation_order_ser');

fprintf('All outputs saved to:\n%s\n', result_dir);

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
    phase = wavenumber * sin(theta) .* ...
        (cos(phi) * y_position + sin(phi) * z_position);
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

function mu = intensity_codebook(G, b, candidates)
    mu = abs(G * candidates.' + b).^2;
end

function [projected_mu, noiseless_power, rank_current, ...
        maximum_relative_residual] = build_reduced_intensity_codebook( ...
        G, b, candidates, cfg, random_seed)
    num_candidates = size(candidates, 1);
    num_observations = size(G, 1);
    mu_reference = intensity_codebook(G, b, candidates(1, :));

    sample_count = min(cfg.basis_sample_size, num_candidates);
    if sample_count == num_candidates
        sample_indices = (1:num_candidates).';
    else
        rng(random_seed, 'twister');
        sample_indices = unique([1; ...
            randperm(num_candidates, sample_count - 1).']);
    end

    mu_sample = intensity_codebook(G, b, candidates(sample_indices, :));
    difference_sample = mu_sample - mu_reference;
    [left_vectors, singular_matrix, ~] = svd(difference_sample, 'econ');
    singular_values = diag(singular_matrix);
    rank_current = sum(singular_values > ...
        cfg.basis_relative_tolerance * max(singular_values));
    basis = left_vectors(:, 1:rank_current);

    projected_mu = zeros(rank_current, num_candidates);
    total_squared_intensity = 0;
    maximum_relative_residual = 0;
    for first_index = 1:cfg.codebook_chunk_size:num_candidates
        last_index = min(first_index + cfg.codebook_chunk_size - 1, ...
            num_candidates);
        indices = first_index:last_index;
        mu_chunk = intensity_codebook(G, b, candidates(indices, :));
        difference_chunk = mu_chunk - mu_reference;
        projected_chunk = basis.' * difference_chunk;
        projected_mu(:, indices) = projected_chunk;
        residual_chunk = difference_chunk - basis * projected_chunk;
        relative_residual = norm(residual_chunk, 'fro') / ...
            max(norm(difference_chunk, 'fro'), eps);
        maximum_relative_residual = max( ...
            maximum_relative_residual, relative_residual);
        total_squared_intensity = total_squared_intensity + ...
            sum(mu_chunk(:).^2);
    end

    noiseless_power = total_squared_intensity / ...
        (num_observations * num_candidates);
    if maximum_relative_residual > cfg.projection_residual_tolerance
        error('Reduced intensity subspace failed the residual check.');
    end
end

function error_count = simulate_exact_projected_ml_errors( ...
        projected_mu, combination_indices, noise_variance, ...
        num_trials, chunk_size, random_seed)
    rng(random_seed, 'twister');
    num_candidates = size(projected_mu, 2);
    reduced_dimension = size(projected_mu, 1);
    num_users = size(combination_indices, 2);
    candidate_norms = sum(projected_mu.^2, 1).';
    error_count = 0;
    completed = 0;

    while completed < num_trials
        current_chunk = min(chunk_size, num_trials - completed);
        true_ids = randi(num_candidates, current_chunk, 1);
        observations = projected_mu(:, true_ids) + ...
            sqrt(noise_variance) * randn(reduced_dimension, current_chunk);
        observation_norms = sum(observations.^2, 1);
        metrics = candidate_norms + observation_norms - ...
            2 * real(projected_mu.' * observations);
        [~, detected_ids] = min(metrics, [], 1);
        detected_ids = detected_ids(:);
        errors = combination_indices(detected_ids, :) ~= ...
            combination_indices(true_ids, :);
        error_count = error_count + sum(errors(:));
        completed = completed + current_chunk;
    end
end

function [lower, upper] = wilson_interval(error_count, sample_count, z)
    p_hat = error_count / sample_count;
    denominator = 1 + z^2 / sample_count;
    center = (p_hat + z^2 / (2 * sample_count)) / denominator;
    half_width = z / denominator * sqrt( ...
        p_hat * (1 - p_hat) / sample_count + ...
        z^2 / (4 * sample_count^2));
    lower = max(center - half_width, 0);
    upper = min(center + half_width, 1);
end

function fig = publication_figure(position)
    fig = figure('Color', 'w', 'Position', position, 'Visible', 'off');
    set(groot, 'defaultAxesFontName', 'Times New Roman', ...
        'defaultTextFontName', 'Times New Roman', ...
        'defaultAxesFontSize', 11, 'defaultAxesLineWidth', 0.8, ...
        'defaultLineLineWidth', 1.5, 'defaultLineMarkerSize', 5.5);
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
