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

%% Configuration inherited from the original perfect-CSI experiment
rng(20260903, 'twister');

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

cfg.num_trials = 20000;
cfg.chunk_size = 2000;

snr_dB_vec = -15:2:3;
reference_names = {'No reference', 'Constant', 'Random', ...
    'Uniform', 'Optimized'};

qpsk = [1 + 1j; -1 + 1j; -1 - 1j; 1 - 1j] / sqrt(2);
H_true = generate_channel_matrix(cfg);
G_data = repmat(H_true, cfg.num_reference_states, 1);
[qpsk_indices, qpsk_candidates] = ...
    enumerate_symbol_vectors(qpsk, cfg.num_users);

[reference_vectors, reference_codebook, optimized_index] = ...
    construct_reference_designs(G_data, qpsk_candidates, cfg);

num_designs = numel(reference_names);
num_snr = numel(snr_dB_vec);
mu_cells = cell(num_designs, 1);
dmin_normalized = zeros(num_designs, 1);

for i_design = 1:num_designs
    mu_cells{i_design} = intensity_codebook( ...
        G_data, reference_vectors(:, i_design), qpsk_candidates);
    dmin = minimum_energy_distance(mu_cells{i_design});
    dmin_normalized(i_design) = dmin / ...
        max(sqrt(mean(mu_cells{i_design}(:).^2)), eps);
end

uniform_equals_optimized = norm( ...
    reference_vectors(:, 4) - reference_vectors(:, 5)) <= 1e-12;

%% Paired high-Monte-Carlo ML detection
% The same transmitted symbol indices and standard-normal noise samples are
% used across reference designs at each SNR. Detection is processed in chunks
% to avoid storing a large observation matrix.
symbol_error_count = zeros(num_snr, num_designs);
ser = zeros(num_snr, num_designs);
num_symbol_decisions = cfg.num_trials * cfg.num_users;

fprintf('\nPerfect-CSI reference-design SER simulation\n');
fprintf('RIS: %d x %d, users: %d, reference states: %d\n', ...
    cfg.ris_rows, cfg.ris_cols, cfg.num_users, ...
    cfg.num_reference_states);
fprintf('Trials per SNR: %d (%d symbol decisions)\n', ...
    cfg.num_trials, num_symbol_decisions);

for i_snr = 1:num_snr
    active_designs = 1:num_designs;
    if uniform_equals_optimized
        active_designs = 1:4;
    end

    error_counts_active = simulate_paired_ml_errors( ...
        mu_cells, active_designs, qpsk_indices, snr_dB_vec(i_snr), ...
        cfg.num_trials, cfg.chunk_size, 610000 + i_snr);
    symbol_error_count(i_snr, active_designs) = error_counts_active;

    if uniform_equals_optimized
        symbol_error_count(i_snr, 5) = symbol_error_count(i_snr, 4);
    end

    ser(i_snr, :) = symbol_error_count(i_snr, :) / ...
        num_symbol_decisions;
    fprintf('SNR %5.1f dB completed\n', snr_dB_vec(i_snr));
end

%% Wilson confidence intervals and plotting floor for zero-error points
ser_ci_lower = zeros(size(ser));
ser_ci_upper = zeros(size(ser));
for i_design = 1:num_designs
    for i_snr = 1:num_snr
        [ser_ci_lower(i_snr, i_design), ...
            ser_ci_upper(i_snr, i_design)] = wilson_interval( ...
            symbol_error_count(i_snr, i_design), ...
            num_symbol_decisions, 1.96);
    end
end

zero_error_floor = 0.5 / num_symbol_decisions;
ser_plot = max(ser, zero_error_floor);
zero_error_indicator = symbol_error_count == 0;

[snr_grid, design_grid] = ndgrid(snr_dB_vec(:), 1:num_designs);
design_name_column = reshape( ...
    string(reference_names(design_grid(:))), [], 1);
summary_table = table( ...
    snr_grid(:), design_name_column, ...
    symbol_error_count(:), ser(:), ser_plot(:), ...
    zero_error_indicator(:), ser_ci_lower(:), ser_ci_upper(:), ...
    dmin_normalized(design_grid(:)), ...
    'VariableNames', {'SNR_dB', 'ReferenceDesign', ...
    'SymbolErrorCount', 'EmpiricalSER', 'PlottedSER', ...
    'ZeroErrorPoint', 'SER_95CI_Lower', 'SER_95CI_Upper', ...
    'NormalizedMinimumEnergyDistance'});
writetable(summary_table, fullfile(result_dir, ...
    'perfect_csi_reference_ser_summary.csv'));

%% Single-panel publication figure
% Match the color, line, and marker mapping used by the channel-NMSE
% reference-design figure in the preceding subsection.
colors = [ ...
    0.20, 0.20, 0.20; ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980; ...
    0.4660, 0.6740, 0.1880];
line_styles = {':', '--', '-.', '-'};
markers = {'x', 'o', 'd', 's'};
plot_names = {'No reference', 'Constant', 'Random', ...
    'Uniform / optimized'};

fig = publication_figure([100, 100, 570, 410]);
hold on;
for i_design = 1:4
    semilogy(snr_dB_vec, ser_plot(:, i_design), ...
        'Color', colors(i_design, :), ...
        'LineStyle', line_styles{i_design}, ...
        'Marker', markers{i_design}, ...
        'LineWidth', 1.5, 'MarkerSize', 5.5, ...
        'MarkerFaceColor', 'w');
end
set(gca, 'YScale', 'log');
grid on;
box on;
xlabel('SNR (dB)');
ylabel('ML SER');
xlim([snr_dB_vec(1), snr_dB_vec(end)]);
ylim([1e-5, 1]);
xticks([-15, -12, -9, -6, -3, 0, 3]);
legend(plot_names, 'Location', 'southwest');
export_publication_figure(fig, result_dir, ...
    'fig_reference_design_ser_perfect_csi');

save(fullfile(result_dir, ...
    'chapter5_perfect_csi_reference_ser_results.mat'), ...
    'cfg', 'snr_dB_vec', 'reference_names', 'reference_vectors', ...
    'reference_codebook', 'optimized_index', ...
    'ser', 'ser_plot', 'symbol_error_count', ...
    'ser_ci_lower', 'ser_ci_upper', 'zero_error_floor', ...
    'zero_error_indicator', 'dmin_normalized', ...
    'uniform_equals_optimized', 'H_true', ...
    'qpsk_indices', 'qpsk_candidates');

fprintf('\nNormalized minimum energy distances:\n');
for i_design = 1:num_designs
    fprintf('%-12s: %.6f\n', ...
        reference_names{i_design}, dmin_normalized(i_design));
end
fprintf('Zero-error plotting floor: %.6e\n', zero_error_floor);
fprintf('Uniform and optimized references identical: %d\n', ...
    uniform_equals_optimized);
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
    phase = wavenumber * (sin(theta) * cos(phi) * y_position + ...
        sin(theta) * sin(phi) * z_position);
    v = exp(1j * phase);
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
        support_fraction = 0.08 + 0.92 * ...
            (i_code - 4) / max(cfg.reference_codebook_size - 4, 1);
        num_active = max(rows, round(support_fraction * num_observations));
        active_indices = randperm(num_observations, num_active);
        active_amplitude = amplitude * ...
            sqrt(num_observations / num_active);
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

function error_counts = simulate_paired_ml_errors( ...
        mu_cells, active_designs, combination_indices, snr_dB, ...
        num_trials, chunk_size, random_seed)
    rng(random_seed, 'twister');
    num_candidates = size(mu_cells{1}, 2);
    num_users = size(combination_indices, 2);
    num_observations = size(mu_cells{1}, 1);
    error_counts = zeros(1, numel(active_designs));
    completed = 0;

    while completed < num_trials
        current_chunk = min(chunk_size, num_trials - completed);
        true_ids = randi(num_candidates, current_chunk, 1);
        standardized_noise = randn(num_observations, current_chunk);

        for i_active = 1:numel(active_designs)
            i_design = active_designs(i_active);
            mu = mu_cells{i_design};
            noise_variance = mean(mu(:).^2) / 10^(snr_dB / 10);
            observations = mu(:, true_ids) + ...
                sqrt(noise_variance) * standardized_noise;

            mu_norms = sum(mu.^2, 1).';
            observation_norms = sum(observations.^2, 1);
            metrics = mu_norms + observation_norms - ...
                2 * real(mu' * observations);
            [~, detected_ids] = min(metrics, [], 1);
            detected_ids = detected_ids(:);
            errors = combination_indices(detected_ids, :) ~= ...
                combination_indices(true_ids, :);
            error_counts(i_active) = error_counts(i_active) + ...
                sum(errors(:));
        end
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
