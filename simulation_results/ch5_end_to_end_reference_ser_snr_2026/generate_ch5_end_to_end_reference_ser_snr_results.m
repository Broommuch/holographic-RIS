%% End-to-end reference-design comparison versus SNR
% This driver is adapted from the Chapter V estimated-CSI simulation. It
% isolates the roles of the training and data references: the training SNR is
% swept at a fixed data SNR in the first panel, whereas the data SNR is swept
% over one common estimated-channel ensemble in the second panel.

clc;
clear;
close all;

script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

rng(20260831, 'twister');

%% Baseline configuration
cfg.num_users = 2;
cfg.ris_rows = 8;
cfg.ris_cols = 8;
cfg.num_ris = cfg.ris_rows * cfg.ris_cols;
cfg.lambda = 1;
cfg.spacing = cfg.lambda / 2;
cfg.theta_true_deg = [15; 35];
cfg.phi_true_deg = [8; -10];
cfg.path_gains = [exp(1j * 0.35); 0.85 * exp(-1j * 0.55)];

cfg.num_pilots = 32;
cfg.reference_amplitude = 1.5;
cfg.training_codebook_size = 32;
cfg.num_reference_states = 4;
cfg.data_codebook_size = 32;

cfg.channel_gn_max_iter = 45;
cfg.channel_gn_tolerance = 1e-7;
cfg.channel_gn_initial_damping = 1e-2;
cfg.channel_gn_num_restarts = 2;
cfg.data_gn_max_iter = 35;
cfg.data_gn_tolerance = 1e-7;
cfg.data_gn_initial_damping = 1e-2;
cfg.data_gs_max_iter = 30;
cfg.data_gs_regularization = 1e-3;

cfg.num_channel_trials_detector = 12;
cfg.num_data_per_channel_detector = 100;
cfg.num_channel_trials_sweep = 15;
cfg.num_data_per_channel_sweep = 120;
cfg.num_controlled_channels = 12;
cfg.num_data_per_controlled_channel = 160;

snr_detector_dB = -12:3:12;
pilot_length_vec = [4, 6, 8, 12, 16, 24, 32, 48];
controlled_nmse_dB = [-30, -25, -20, -16, -13, -10, -7, -4, -1, 2];
reference_names = {'No reference', 'Constant', 'Random', ...
    'Uniform', 'Optimized'};

training_reference_snr_vec = -6:3:18;
data_reference_snr_vec = -12:3:12;
cfg.reference_curve_data_snr_dB = -6;
cfg.reference_curve_training_snr_dB = 12;
cfg.num_channel_trials_reference_curve = 20;
cfg.num_data_per_channel_reference_curve = 250;

cfg.pilot_sweep_training_snr_dB = 12;
cfg.pilot_sweep_data_snr_dB = -6;
cfg.controlled_nmse_data_snr_dB = -3;
cfg.training_reference_snr_dB = 6;
cfg.training_reference_data_snr_dB = -6;
cfg.data_reference_training_snr_dB = 12;
cfg.data_reference_snr_dB = -6;

qpsk = [1 + 1j; -1 + 1j; -1 - 1j; 1 - 1j] / sqrt(2);
[combination_indices, symbol_candidates] = ...
    enumerate_symbol_vectors(qpsk, cfg.num_users);
H_true = generate_channel_matrix(cfg);
G_true = repmat(H_true, cfg.num_reference_states, 1);

max_pilots = max([cfg.num_pilots, pilot_length_vec]);
pilot_symbols_full = generate_qpsk_pilots( ...
    max_pilots, cfg.num_users, 20260901);
pilot_symbols = pilot_symbols_full(1:cfg.num_pilots, :);
training_phase_optimized = select_optimized_training_reference_phase( ...
    pilot_symbols, H_true, cfg.training_codebook_size, 20260902);
B_training_baseline = cfg.reference_amplitude * training_phase_optimized;

[data_reference_vectors, data_reference_dmin, data_optimized_index] = ...
    construct_data_reference_designs(G_true, symbol_candidates, cfg);
b_data_baseline = data_reference_vectors(:, 5);

fprintf('\nChapter V-D end-to-end detection with estimated CSI\n');
fprintf('RIS: %d x %d, users: %d, pilots: %d\n', ...
    cfg.ris_rows, cfg.ris_cols, cfg.num_users, cfg.num_pilots);
fprintf('Selected data-reference codeword: %d of %d\n', ...
    data_optimized_index, cfg.data_codebook_size);

%% 1) Perfect versus estimated CSI for three detectors
if false
num_snr = numel(snr_detector_dB);
ser_perfect_ml = zeros(num_snr, 1);
ser_perfect_gn = zeros(num_snr, 1);
ser_perfect_gs = zeros(num_snr, 1);
ser_estimated_ml = zeros(num_snr, 1);
ser_estimated_gn = zeros(num_snr, 1);
ser_estimated_gs = zeros(num_snr, 1);
channel_nmse_detector_dB = zeros(num_snr, 1);

num_perfect_trials = cfg.num_channel_trials_detector * ...
    cfg.num_data_per_channel_detector;

for i_snr = 1:num_snr
    snr_dB = snr_detector_dB(i_snr);
    [ser_perfect_ml(i_snr), ser_perfect_gn(i_snr), ...
        ser_perfect_gs(i_snr)] = simulate_three_detector_ser( ...
            G_true, G_true, b_data_baseline, qpsk, ...
            combination_indices, symbol_candidates, snr_dB, ...
            num_perfect_trials, 10000 + i_snr, cfg);

    ml_trials = zeros(cfg.num_channel_trials_detector, 1);
    gn_trials = zeros(cfg.num_channel_trials_detector, 1);
    gs_trials = zeros(cfg.num_channel_trials_detector, 1);
    nmse_trials = zeros(cfg.num_channel_trials_detector, 1);

    for channel_trial = 1:cfg.num_channel_trials_detector
        H_hat = estimate_channel_realization( ...
            pilot_symbols, B_training_baseline, H_true, snr_dB, ...
            11000 + 100 * i_snr + channel_trial, cfg);
        nmse_trials(channel_trial) = channel_nmse(H_hat, H_true);
        G_hat = repmat(H_hat, cfg.num_reference_states, 1);
        [ml_trials(channel_trial), gn_trials(channel_trial), ...
            gs_trials(channel_trial)] = simulate_three_detector_ser( ...
                G_true, G_hat, b_data_baseline, qpsk, ...
                combination_indices, symbol_candidates, snr_dB, ...
                cfg.num_data_per_channel_detector, ...
                12000 + 1000 * i_snr + channel_trial, cfg);
    end

    ser_estimated_ml(i_snr) = mean(ml_trials);
    ser_estimated_gn(i_snr) = mean(gn_trials);
    ser_estimated_gs(i_snr) = mean(gs_trials);
    channel_nmse_detector_dB(i_snr) = ...
        10 * log10(max(mean(nmse_trials), eps));

    fprintf(['SNR %5.1f dB: perfect ML/GN/GS %.3e/%.3e/%.3e, ', ...
        'estimated %.3e/%.3e/%.3e, NMSE %.2f dB\n'], ...
        snr_dB, ser_perfect_ml(i_snr), ser_perfect_gn(i_snr), ...
        ser_perfect_gs(i_snr), ser_estimated_ml(i_snr), ...
        ser_estimated_gn(i_snr), ser_estimated_gs(i_snr), ...
        channel_nmse_detector_dB(i_snr));
end

detector_table = table(snr_detector_dB(:), channel_nmse_detector_dB, ...
    ser_perfect_ml, ser_estimated_ml, ser_perfect_gn, ser_estimated_gn, ...
    ser_perfect_gs, ser_estimated_gs, ...
    'VariableNames', {'SNR_dB', 'Channel_NMSE_dB', ...
    'Perfect_ML_SER', 'Estimated_ML_SER', ...
    'Perfect_GN_SER', 'Estimated_GN_SER', ...
    'Perfect_GS_SER', 'Estimated_GS_SER'});
writetable(detector_table, ...
    fullfile(result_dir, 'perfect_estimated_csi_detector_ser.csv'));

fig = publication_figure([100, 100, 570, 410]);
colors = lines(3);
semilogy(snr_detector_dB, plotting_ser(ser_perfect_ml, ...
    num_perfect_trials, cfg.num_users), 'o--', 'Color', colors(1, :));
hold on;
semilogy(snr_detector_dB, plotting_ser(ser_estimated_ml, ...
    num_perfect_trials, cfg.num_users), 'o-', 'Color', colors(1, :), ...
    'MarkerFaceColor', colors(1, :));
semilogy(snr_detector_dB, plotting_ser(ser_perfect_gn, ...
    num_perfect_trials, cfg.num_users), 's--', 'Color', colors(2, :));
semilogy(snr_detector_dB, plotting_ser(ser_estimated_gn, ...
    num_perfect_trials, cfg.num_users), 's-', 'Color', colors(2, :), ...
    'MarkerFaceColor', colors(2, :));
semilogy(snr_detector_dB, plotting_ser(ser_perfect_gs, ...
    num_perfect_trials, cfg.num_users), 'd--', 'Color', colors(3, :));
semilogy(snr_detector_dB, plotting_ser(ser_estimated_gs, ...
    num_perfect_trials, cfg.num_users), 'd-', 'Color', colors(3, :), ...
    'MarkerFaceColor', colors(3, :));
grid on;
box on;
xlabel('Training/data SNR (dB)');
ylabel('SER');
legend('ML, perfect CSI', 'ML, estimated CSI', ...
    'GN, perfect CSI', 'GN, estimated CSI', ...
    'GS, perfect CSI', 'GS, estimated CSI', ...
    'Location', 'southwest', 'NumColumns', 2);
xlim([snr_detector_dB(1), snr_detector_dB(end)]);
ylim([2e-4, 1]);
export_publication_figure(fig, result_dir, ...
    'fig_ser_perfect_estimated_csi');

%% 2) Pilot length and controlled channel mismatch
num_pilot_lengths = numel(pilot_length_vec);
ser_pilot_length = zeros(num_pilot_lengths, 1);
channel_nmse_pilot_dB = zeros(num_pilot_lengths, 1);

for i_length = 1:num_pilot_lengths
    num_pilots = pilot_length_vec(i_length);
    S = pilot_symbols_full(1:num_pilots, :);
    phase_reference = select_optimized_training_reference_phase( ...
        S, H_true, max(12, floor(cfg.training_codebook_size / 2)), ...
        20000 + i_length);
    B = cfg.reference_amplitude * phase_reference;
    ser_trials = zeros(cfg.num_channel_trials_sweep, 1);
    nmse_trials = zeros(cfg.num_channel_trials_sweep, 1);

    for channel_trial = 1:cfg.num_channel_trials_sweep
        H_hat = estimate_channel_realization( ...
            S, B, H_true, cfg.pilot_sweep_training_snr_dB, ...
            21000 + 100 * i_length + channel_trial, cfg);
        nmse_trials(channel_trial) = channel_nmse(H_hat, H_true);
        G_hat = repmat(H_hat, cfg.num_reference_states, 1);
        ser_trials(channel_trial) = simulate_ml_mismatched_ser( ...
            G_true, G_hat, b_data_baseline, combination_indices, ...
            symbol_candidates, cfg.pilot_sweep_data_snr_dB, ...
            cfg.num_data_per_channel_sweep, ...
            22000 + 1000 * i_length + channel_trial);
    end

    ser_pilot_length(i_length) = mean(ser_trials);
    channel_nmse_pilot_dB(i_length) = ...
        10 * log10(max(mean(nmse_trials), eps));
    fprintf('Pilots %2d: NMSE %.2f dB, end-to-end ML SER %.3e\n', ...
        num_pilots, channel_nmse_pilot_dB(i_length), ...
        ser_pilot_length(i_length));
end

pilot_table = table(pilot_length_vec(:), channel_nmse_pilot_dB, ...
    ser_pilot_length, 'VariableNames', ...
    {'PilotLength', 'Channel_NMSE_dB', 'EndToEnd_ML_SER'});
writetable(pilot_table, fullfile(result_dir, 'pilot_length_end_to_end.csv'));

num_nmse = numel(controlled_nmse_dB);
ser_controlled_nmse = zeros(num_nmse, 1);
actual_controlled_nmse_dB = zeros(num_nmse, 1);

for i_nmse = 1:num_nmse
    ser_trials = zeros(cfg.num_controlled_channels, 1);
    actual_nmse_trials = zeros(cfg.num_controlled_channels, 1);
    for channel_trial = 1:cfg.num_controlled_channels
        H_hat = perturb_channel_to_nmse( ...
            H_true, controlled_nmse_dB(i_nmse), ...
            30000 + 100 * i_nmse + channel_trial);
        actual_nmse_trials(channel_trial) = channel_nmse(H_hat, H_true);
        G_hat = repmat(H_hat, cfg.num_reference_states, 1);
        ser_trials(channel_trial) = simulate_ml_mismatched_ser( ...
            G_true, G_hat, b_data_baseline, combination_indices, ...
            symbol_candidates, cfg.controlled_nmse_data_snr_dB, ...
            cfg.num_data_per_controlled_channel, ...
            31000 + 1000 * i_nmse + channel_trial);
    end
    ser_controlled_nmse(i_nmse) = mean(ser_trials);
    actual_controlled_nmse_dB(i_nmse) = ...
        10 * log10(mean(actual_nmse_trials));
end

controlled_table = table(actual_controlled_nmse_dB, ser_controlled_nmse, ...
    'VariableNames', {'Channel_NMSE_dB', 'EndToEnd_ML_SER'});
writetable(controlled_table, fullfile(result_dir, 'controlled_nmse_ser.csv'));

fig = publication_figure([100, 100, 560, 650]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
yyaxis left;
semilogy(pilot_length_vec, plotting_ser(ser_pilot_length, ...
    cfg.num_channel_trials_sweep * cfg.num_data_per_channel_sweep, ...
    cfg.num_users), 'o-', 'Color', colors(1, :));
ylabel('End-to-end ML SER');
ylim([2e-4, 1]);
yyaxis right;
plot(pilot_length_vec, channel_nmse_pilot_dB, 's--', ...
    'Color', colors(2, :));
ylabel('Channel NMSE (dB)');
xlabel('Pilot length T_p');
grid on;
box on;
title('(a) Pilot-length tradeoff');

nexttile;
semilogy(actual_controlled_nmse_dB, plotting_ser(ser_controlled_nmse, ...
    cfg.num_controlled_channels * cfg.num_data_per_controlled_channel, ...
    cfg.num_users), 'o-', 'Color', colors(1, :), ...
    'MarkerFaceColor', colors(1, :));
grid on;
box on;
xlabel('Channel NMSE (dB)');
ylabel('End-to-end ML SER');
ylim([2e-4, 1]);
title('(b) Controlled CSI mismatch');
export_publication_figure(fig, result_dir, ...
    'fig_pilot_nmse_end_to_end');

end

%% 3) Training- and data-phase reference designs
num_designs = numel(reference_names);
training_reference_cells = construct_training_reference_designs( ...
    pilot_symbols, H_true, cfg);
num_training_snr = numel(training_reference_snr_vec);
num_data_snr = numel(data_reference_snr_vec);
num_channel_trials = cfg.num_channel_trials_reference_curve;
num_data_trials = cfg.num_data_per_channel_reference_curve;
num_symbol_decisions = num_channel_trials * num_data_trials * ...
    cfg.num_users;
zero_error_floor = 0.5 / num_symbol_decisions;

training_error_count = zeros(num_training_snr, num_designs);
training_ser = zeros(num_training_snr, num_designs);
training_nmse_dB = zeros(num_training_snr, num_designs);

fprintf('\nTraining-reference SER curves at data SNR %.1f dB\n', ...
    cfg.reference_curve_data_snr_dB);
for i_snr = 1:num_training_snr
    training_snr_dB = training_reference_snr_vec(i_snr);
    for i_design = 1:num_designs
        B = training_reference_cells{i_design};
        nmse_trials = zeros(num_channel_trials, 1);
        for channel_trial = 1:num_channel_trials
            % Reusing the seed across designs and SNRs pairs the pilot-noise
            % realizations without changing any marginal distribution.
            H_hat = estimate_channel_realization( ...
                pilot_symbols, B, H_true, training_snr_dB, ...
                610000 + channel_trial, cfg);
            nmse_trials(channel_trial) = channel_nmse(H_hat, H_true);
            G_hat = repmat(H_hat, cfg.num_reference_states, 1);
            ser_trial = simulate_ml_mismatched_ser( ...
                G_true, G_hat, b_data_baseline, combination_indices, ...
                symbol_candidates, cfg.reference_curve_data_snr_dB, ...
                num_data_trials, 620000 + channel_trial);
            training_error_count(i_snr, i_design) = ...
                training_error_count(i_snr, i_design) + round( ...
                ser_trial * num_data_trials * cfg.num_users);
        end
        training_ser(i_snr, i_design) = ...
            training_error_count(i_snr, i_design) / num_symbol_decisions;
        training_nmse_dB(i_snr, i_design) = 10 * log10( ...
            max(mean(nmse_trials), eps));
    end
    fprintf('Training SNR %5.1f dB completed\n', training_snr_dB);
end

training_ser_plot = max(training_ser, zero_error_floor);
[training_snr_grid, training_design_grid] = ndgrid( ...
    training_reference_snr_vec(:), 1:num_designs);
training_design_column = reshape( ...
    string(reference_names(training_design_grid(:))), [], 1);
training_reference_table = table( ...
    training_snr_grid(:), training_design_column, ...
    training_error_count(:), training_ser(:), training_ser_plot(:), ...
    training_nmse_dB(:), ...
    'VariableNames', {'TrainingSNR_dB', 'TrainingReference', ...
    'SymbolErrorCount', 'EmpiricalSER', 'PlottedSER', ...
    'ChannelNMSE_dB'});
writetable(training_reference_table, fullfile(result_dir, ...
    'training_reference_ser_snr.csv'));

% Estimate one common channel ensemble, then change only the data reference
% and data SNR. This prevents channel-estimation variability from obscuring
% the comparison among data-reference designs.
common_channel_estimates = cell(num_channel_trials, 1);
for channel_trial = 1:num_channel_trials
    common_channel_estimates{channel_trial} = estimate_channel_realization( ...
        pilot_symbols, B_training_baseline, H_true, ...
        cfg.reference_curve_training_snr_dB, ...
        630000 + channel_trial, cfg);
end

data_error_count = zeros(num_data_snr, num_designs);
data_ser = zeros(num_data_snr, num_designs);
fprintf('\nData-reference SER curves at training SNR %.1f dB\n', ...
    cfg.reference_curve_training_snr_dB);
for i_snr = 1:num_data_snr
    data_snr_dB = data_reference_snr_vec(i_snr);
    for i_design = 1:num_designs
        for channel_trial = 1:num_channel_trials
            H_hat = common_channel_estimates{channel_trial};
            G_hat = repmat(H_hat, cfg.num_reference_states, 1);
            ser_trial = simulate_ml_mismatched_ser( ...
                G_true, G_hat, data_reference_vectors(:, i_design), ...
                combination_indices, symbol_candidates, data_snr_dB, ...
                num_data_trials, 640000 + channel_trial);
            data_error_count(i_snr, i_design) = ...
                data_error_count(i_snr, i_design) + round( ...
                ser_trial * num_data_trials * cfg.num_users);
        end
        data_ser(i_snr, i_design) = ...
            data_error_count(i_snr, i_design) / num_symbol_decisions;
    end
    fprintf('Data SNR %5.1f dB completed\n', data_snr_dB);
end

uniform_equals_optimized = norm( ...
    data_reference_vectors(:, 4) - data_reference_vectors(:, 5)) <= 1e-12;
data_ser_plot = max(data_ser, zero_error_floor);
[data_snr_grid, data_design_grid] = ndgrid( ...
    data_reference_snr_vec(:), 1:num_designs);
data_design_column = reshape( ...
    string(reference_names(data_design_grid(:))), [], 1);
data_reference_table = table( ...
    data_snr_grid(:), data_design_column, data_error_count(:), ...
    data_ser(:), data_ser_plot(:), ...
    data_reference_dmin(data_design_grid(:)), ...
    'VariableNames', {'DataSNR_dB', 'DataReference', ...
    'SymbolErrorCount', 'EmpiricalSER', 'PlottedSER', ...
    'NormalizedMinimumEnergyDistance'});
writetable(data_reference_table, fullfile(result_dir, ...
    'data_reference_ser_snr.csv'));

% The first five styles match the mapping used in Figs. 5 and 7. The
% optimized training reference needs the fifth (purple triangle) style,
% whereas the uniform and optimized data references are identical and are
% therefore represented by one green curve in the second panel.
colors = [ ...
    0.20, 0.20, 0.20; ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980; ...
    0.4660, 0.6740, 0.1880; ...
    0.4940, 0.1840, 0.5560];
line_styles = {':', '--', '-.', '-', '-'};
markers = {'x', 'o', 'd', 's', '^'};

fig = publication_figure([100, 100, 590, 660]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
hold on;
for i_design = 1:num_designs
    semilogy(training_reference_snr_vec, ...
        training_ser_plot(:, i_design), ...
        'Color', colors(i_design, :), ...
        'LineStyle', line_styles{i_design}, ...
        'Marker', markers{i_design}, 'LineWidth', 1.5, ...
        'MarkerSize', 5.5, 'MarkerFaceColor', 'w');
end
set(gca, 'YScale', 'log');
grid on;
box on;
ylabel('End-to-end ML SER');
xlabel('Training SNR (dB)');
xlim([training_reference_snr_vec(1), training_reference_snr_vec(end)]);
ylim([3e-5, 1]);
legend(reference_names, 'Location', 'southwest', 'NumColumns', 2);
title('(a) Training-reference design');

nexttile;
hold on;
for i_design = 1:4
    semilogy(data_reference_snr_vec, data_ser_plot(:, i_design), ...
        'Color', colors(i_design, :), ...
        'LineStyle', line_styles{i_design}, ...
        'Marker', markers{i_design}, 'LineWidth', 1.5, ...
        'MarkerSize', 5.5, 'MarkerFaceColor', 'w');
end
set(gca, 'YScale', 'log');
grid on;
box on;
ylabel('End-to-end ML SER');
xlabel('Data SNR (dB)');
xlim([data_reference_snr_vec(1), data_reference_snr_vec(end)]);
ylim([3e-5, 1]);
legend({'No reference', 'Constant', 'Random', ...
    'Uniform / optimized'}, 'Location', 'southwest');
title('(b) Data-reference design');
export_publication_figure(fig, result_dir, ...
    'fig_reference_design_ser_snr_end_to_end');

%% Save complete numerical workspace
save(fullfile(result_dir, ...
    'chapter5_end_to_end_reference_ser_snr_results.mat'), ...
    'cfg', 'training_reference_snr_vec', 'data_reference_snr_vec', ...
    'reference_names', 'training_reference_cells', ...
    'training_error_count', 'training_ser', 'training_ser_plot', ...
    'training_nmse_dB', 'data_reference_vectors', ...
    'data_reference_dmin', 'data_error_count', 'data_ser', ...
    'data_ser_plot', 'zero_error_floor', ...
    'uniform_equals_optimized', 'H_true', 'B_training_baseline');

fprintf('\nAll outputs saved to:\n%s\n', result_dir);

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
    for user = 1:cfg.num_users
        H(:, user) = cfg.path_gains(user) * steering_vector( ...
            cfg.theta_true_deg(user), cfg.phi_true_deg(user), cfg);
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

function phase_reference = make_random_training_reference_phase( ...
        num_pilots, num_ris, seed)
    old_rng = rng;
    rng(seed, 'twister');
    phase_reference = exp(1j * 2 * pi * rand(num_pilots, num_ris));
    rng(old_rng);
end

function phase_reference = make_uniform_training_reference_phase( ...
        num_pilots, num_ris)
    time_phase = 2 * pi * (0:num_pilots - 1).' / num_pilots;
    space_phase = 2 * pi * (0:num_ris - 1) / num_ris;
    phase_reference = exp(1j * (time_phase + space_phase));
end

function phase_best = select_optimized_training_reference_phase( ...
        S, H, num_candidates, seed)
    num_pilots = size(S, 1);
    num_ris = size(H, 1);
    phase_best = make_uniform_training_reference_phase(num_pilots, num_ris);
    best_score = channel_jacobian_minimum_singular_value(S, H, phase_best);
    for candidate = 1:num_candidates
        phase_candidate = make_random_training_reference_phase( ...
            num_pilots, num_ris, seed + candidate);
        score = channel_jacobian_minimum_singular_value(S, H, phase_candidate);
        if score > best_score
            best_score = score;
            phase_best = phase_candidate;
        end
    end
end

function sigma_min_global = channel_jacobian_minimum_singular_value(S, H, B)
    sigma_min_global = inf;
    for element = 1:size(H, 1)
        h = H(element, :).';
        q = S * h + B(:, element);
        weighted = conj(q) .* S;
        jacobian = 2 * [real(weighted), -imag(weighted)];
        singular_values = svd(jacobian, 'econ');
        sigma_min_global = min(sigma_min_global, min(singular_values));
    end
end

function references = construct_training_reference_designs(S, H, cfg)
    num_pilots = size(S, 1);
    references = cell(5, 1);
    references{1} = zeros(num_pilots, cfg.num_ris);
    references{2} = cfg.reference_amplitude * ...
        ones(num_pilots, cfg.num_ris);
    references{3} = cfg.reference_amplitude * ...
        make_random_training_reference_phase( ...
        num_pilots, cfg.num_ris, 20260903);
    references{4} = cfg.reference_amplitude * ...
        make_uniform_training_reference_phase(num_pilots, cfg.num_ris);
    references{5} = cfg.reference_amplitude * ...
        select_optimized_training_reference_phase( ...
        S, H, cfg.training_codebook_size, 20260904);
end

function [reference_vectors, normalized_dmin, optimized_index] = ...
        construct_data_reference_designs(G, candidates, cfg)
    states = cfg.num_reference_states;
    elements = cfg.num_ris;
    amplitude = cfg.reference_amplitude;
    b_none = zeros(states * elements, 1);
    b_constant = amplitude * ones(states * elements, 1);

    old_rng = rng;
    rng(20260905, 'twister');
    random_matrix = amplitude * exp(1j * 2 * pi * rand(states, elements));
    b_random = reshape(random_matrix.', [], 1);
    state_phase = (0:states - 1).' * 2 * pi / states;
    uniform_matrix = amplitude * exp(1j * state_phase) * ones(1, elements);
    b_uniform = reshape(uniform_matrix.', [], 1);

    codebook = zeros(states * elements, cfg.data_codebook_size);
    codebook(:, 1) = b_uniform;
    codebook(:, 2) = b_constant;
    codebook(:, 3) = b_random;
    for candidate = 4:cfg.data_codebook_size
        phase_matrix = 2 * pi * rand(states, elements);
        codebook(:, candidate) = reshape( ...
            (amplitude * exp(1j * phase_matrix)).', [], 1);
    end
    rng(old_rng);

    codebook_score = zeros(cfg.data_codebook_size, 1);
    for candidate = 1:cfg.data_codebook_size
        mu = abs(G * candidates.' + codebook(:, candidate)).^2;
        codebook_score(candidate) = minimum_energy_distance(mu) / ...
            max(sqrt(mean(mu(:).^2)), eps);
    end
    [~, optimized_index] = max(codebook_score);
    b_optimized = codebook(:, optimized_index);
    reference_vectors = [b_none, b_constant, b_random, ...
        b_uniform, b_optimized];

    normalized_dmin = zeros(5, 1);
    for design = 1:5
        mu = abs(G * candidates.' + reference_vectors(:, design)).^2;
        normalized_dmin(design) = minimum_energy_distance(mu) / ...
            max(sqrt(mean(mu(:).^2)), eps);
    end
end

function distance = minimum_energy_distance(mu)
    squared_norms = sum(mu.^2, 1);
    distance_squared = squared_norms.' + squared_norms - ...
        2 * real(mu' * mu);
    distance_squared = max(distance_squared, 0);
    distance_squared(1:size(distance_squared, 1) + 1:end) = inf;
    distance = sqrt(min(distance_squared(:)));
end

function H_hat = estimate_channel_realization(S, B, H_true, snr_dB, ...
        random_seed, cfg)
    rng(random_seed, 'twister');
    mu = abs(S * H_true.' + B).^2;
    noise_variance = mean(mu(:).^2) / 10^(snr_dB / 10);
    z = mu + sqrt(noise_variance) * randn(size(mu));
    H_hat = recover_multiuser_channel_intensity_gn(z, S, B, cfg);
end

function H_hat = recover_multiuser_channel_intensity_gn(z, S, B, cfg)
    H_hat = zeros(size(z, 2), cfg.num_users);
    for element = 1:size(z, 2)
        H_hat(element, :) = recover_single_element_channel_gn( ...
            z(:, element), S, B(:, element), cfg).';
    end
end

function h_best = recover_single_element_channel_gn(z, S, b, cfg)
    users = size(S, 2);
    z = z(:);
    b = b(:);
    weighted = conj(b) .* S;
    linear_matrix = 2 * [real(weighted), -imag(weighted)];
    linear_target = z - abs(b).^2;
    regularization = 1e-3 * trace(linear_matrix.' * linear_matrix) / ...
        max(2 * users, 1);

    if norm(linear_matrix, 'fro') > 1e-12
        x_linear = (linear_matrix.' * linear_matrix + ...
            (regularization + 1e-9) * eye(2 * users)) \ ...
            (linear_matrix.' * linear_target);
        h_linear = x_linear(1:users) + 1j * x_linear(users + 1:end);
    else
        h_linear = zeros(users, 1);
    end

    signal_scale = sqrt(max(mean(max(z, 0)), 1e-3) / max(users, 1));
    h_best = h_linear;
    best_cost = inf;

    for restart = 1:cfg.channel_gn_num_restarts
        if restart == 1 && norm(h_linear) > 1e-10
            h = h_linear;
        else
            h = signal_scale * (randn(users, 1) + ...
                1j * randn(users, 1)) / sqrt(2);
        end
        damping = cfg.channel_gn_initial_damping;
        current_cost = channel_intensity_cost(h, S, b, z);

        for iteration = 1:cfg.channel_gn_max_iter
            q = S * h + b;
            residual = abs(q).^2 - z;
            weighted = conj(q) .* S;
            jacobian = 2 * [real(weighted), -imag(weighted)];
            normal_matrix = jacobian.' * jacobian;
            gradient = jacobian.' * residual;
            diagonal_scaling = diag(max(diag(normal_matrix), 1e-9));
            accepted = false;
            increment = zeros(2 * users, 1);

            for damping_trial = 1:10
                system_matrix = normal_matrix + damping * diagonal_scaling + ...
                    1e-10 * eye(2 * users);
                trial_increment = -system_matrix \ gradient;
                h_trial = h + trial_increment(1:users) + ...
                    1j * trial_increment(users + 1:end);
                trial_cost = channel_intensity_cost(h_trial, S, b, z);
                if trial_cost < current_cost
                    h = h_trial;
                    current_cost = trial_cost;
                    increment = trial_increment;
                    damping = max(damping / 3, 1e-10);
                    accepted = true;
                    break;
                end
                damping = min(damping * 10, 1e10);
            end

            if ~accepted
                break;
            end
            if norm(increment) <= cfg.channel_gn_tolerance * ...
                    (norm([real(h); imag(h)]) + 1)
                break;
            end
        end

        if current_cost < best_cost
            best_cost = current_cost;
            h_best = h;
        end
    end
end

function value = channel_intensity_cost(h, S, b, z)
    value = mean((abs(S * h + b).^2 - z).^2);
end

function value = channel_nmse(H_hat, H_true)
    value = norm(H_hat - H_true, 'fro')^2 / max(norm(H_true, 'fro')^2, eps);
end

function H_hat = perturb_channel_to_nmse(H_true, target_nmse_dB, seed)
    rng(seed, 'twister');
    perturbation = randn(size(H_true)) + 1j * randn(size(H_true));
    perturbation = perturbation / max(norm(perturbation, 'fro'), eps);
    perturbation = perturbation * norm(H_true, 'fro') * ...
        sqrt(10^(target_nmse_dB / 10));
    H_hat = H_true + perturbation;
end

function ser = simulate_ml_mismatched_ser(G_true, G_hat, b, ...
        combination_indices, candidates, snr_dB, num_trials, seed)
    rng(seed, 'twister');
    mu_true = abs(G_true * candidates.' + b).^2;
    mu_hat = abs(G_hat * candidates.' + b).^2;
    num_candidates = size(candidates, 1);
    users = size(combination_indices, 2);
    true_ids = randi(num_candidates, num_trials, 1);
    noise_variance = mean(mu_true(:).^2) / 10^(snr_dB / 10);
    z = mu_true(:, true_ids) + sqrt(noise_variance) * ...
        randn(size(mu_true, 1), num_trials);
    metrics = sum(mu_hat.^2, 1).' + sum(z.^2, 1) - ...
        2 * real(mu_hat' * z);
    [~, detected_ids] = min(metrics, [], 1);
    detected_ids = detected_ids(:);
    errors = combination_indices(detected_ids, :) ~= ...
        combination_indices(true_ids, :);
    ser = sum(errors(:)) / (num_trials * users);
end

function [ser_ml, ser_gn, ser_gs] = simulate_three_detector_ser( ...
        G_true, G_hat, b, constellation, combination_indices, candidates, ...
        snr_dB, num_trials, seed, cfg)
    rng(seed, 'twister');
    mu_true = abs(G_true * candidates.' + b).^2;
    mu_hat = abs(G_hat * candidates.' + b).^2;
    num_candidates = size(candidates, 1);
    users = size(combination_indices, 2);
    true_ids = randi(num_candidates, num_trials, 1);
    noise_variance = mean(mu_true(:).^2) / 10^(snr_dB / 10);
    z_all = mu_true(:, true_ids) + sqrt(noise_variance) * ...
        randn(size(mu_true, 1), num_trials);

    metrics = sum(mu_hat.^2, 1).' + sum(z_all.^2, 1) - ...
        2 * real(mu_hat' * z_all);
    [~, detected_ml_ids] = min(metrics, [], 1);
    detected_ml = combination_indices(detected_ml_ids(:), :);
    true_indices = combination_indices(true_ids, :);

    detected_gn = zeros(num_trials, users);
    detected_gs = zeros(num_trials, users);
    scale = trace(G_hat' * G_hat) / max(users, 1);
    back_projection = (G_hat' * G_hat + ...
        cfg.data_gs_regularization * max(scale, 1e-6) * eye(users)) \ G_hat';
    for trial = 1:num_trials
        z = z_all(:, trial);
        initial = data_linearized_initialization(z, G_hat, b);
        detected_gn(trial, :) = projected_data_gn( ...
            z, G_hat, b, constellation, initial, cfg).';
        detected_gs(trial, :) = reference_assisted_data_gs( ...
            z, G_hat, b, constellation, initial, back_projection, cfg).';
    end

    ser_ml = sum(detected_ml(:) ~= true_indices(:)) / (num_trials * users);
    ser_gn = sum(detected_gn(:) ~= true_indices(:)) / (num_trials * users);
    ser_gs = sum(detected_gs(:) ~= true_indices(:)) / (num_trials * users);
end

function initial = data_linearized_initialization(z, G, b)
    users = size(G, 2);
    weighted = conj(b) .* G;
    linear_matrix = 2 * [real(weighted), -imag(weighted)];
    target = z - abs(b).^2;
    regularization = 1e-3 * trace(linear_matrix.' * linear_matrix) / ...
        max(2 * users, 1);
    x = (linear_matrix.' * linear_matrix + ...
        (regularization + 1e-9) * eye(2 * users)) \ ...
        (linear_matrix.' * target);
    initial = x(1:users) + 1j * x(users + 1:end);
end

function indices = projected_data_gn(z, G, b, constellation, initial, cfg)
    users = size(G, 2);
    starting_points = [initial, zeros(users, 1)];
    best_cost = inf;
    best_indices = ones(users, 1);
    for restart = 1:size(starting_points, 2)
        s = starting_points(:, restart);
        damping = cfg.data_gn_initial_damping;
        cost = norm(z - abs(G * s + b).^2)^2;
        for iteration = 1:cfg.data_gn_max_iter
            q = G * s + b;
            predicted = abs(q).^2;
            weighted = conj(q) .* G;
            jacobian = 2 * [real(weighted), -imag(weighted)];
            increment = (jacobian.' * jacobian + ...
                damping * eye(2 * users)) \ (jacobian.' * (z - predicted));
            s_trial = s + increment(1:users) + ...
                1j * increment(users + 1:end);
            trial_cost = norm(z - abs(G * s_trial + b).^2)^2;
            if trial_cost < cost
                s = s_trial;
                cost = trial_cost;
                damping = max(damping / 2, 1e-9);
            else
                damping = min(damping * 5, 1e9);
            end
            if norm(increment) / max(norm([real(s); imag(s)]), 1) < ...
                    cfg.data_gn_tolerance
                break;
            end
        end
        trial_indices = quantize_indices(s, constellation);
        projected = constellation(trial_indices);
        projected_cost = norm(z - abs(G * projected + b).^2)^2;
        if projected_cost < best_cost
            best_cost = projected_cost;
            best_indices = trial_indices;
        end
    end
    indices = best_indices;
end

function indices = reference_assisted_data_gs( ...
        z, G, b, constellation, initial, back_projection, cfg)
    users = size(G, 2);
    initial_indices = quantize_indices(initial, constellation);
    starting_points = [constellation(initial_indices), zeros(users, 1)];
    measured_magnitude = sqrt(max(z, 0));
    best_cost = inf;
    best_indices = initial_indices;
    for restart = 1:size(starting_points, 2)
        s = starting_points(:, restart);
        for iteration = 1:cfg.data_gs_max_iter
            q = G * s + b;
            phase = q ./ max(abs(q), eps);
            continuous = back_projection * (measured_magnitude .* phase - b);
            trial_indices = quantize_indices(continuous, constellation);
            s_next = constellation(trial_indices);
            if norm(s_next - s) / max(norm(s), 1) < cfg.data_gn_tolerance
                s = s_next;
                break;
            end
            s = s_next;
        end
        trial_indices = quantize_indices(s, constellation);
        projected = constellation(trial_indices);
        cost = norm(z - abs(G * projected + b).^2)^2;
        if cost < best_cost
            best_cost = cost;
            best_indices = trial_indices;
        end
    end
    indices = best_indices;
end

function indices = quantize_indices(symbols, constellation)
    distances = abs(symbols(:) - constellation(:).');
    [~, indices] = min(distances, [], 2);
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
