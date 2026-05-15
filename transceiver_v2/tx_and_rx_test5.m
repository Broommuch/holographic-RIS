%% 这个脚本准备尝试对多用户进行doa估计

clc; clear; close all;

%% ================= 多用户 DOA 导频发送端参数设置 =================

U = 3;                       % 用户数
N_pilot = 20;               % 每个用户导频符号数
M_mod = 4;                   % QPSK
bits_per_sym = log2(M_mod);

Rs = 1e6;                    % 符号率
sps = 8;                     % 每符号基带采样点
rolloff = 0.25;              % RRC 滚降系数
span = 6;                    % RRC 滤波器跨度，单位：符号

interp = 20;                 % 额外插值倍数
fc = 20e6;                   % 仿真载波频率

fs_bb = Rs * sps;            % 基带采样率
fs = fs_bb * interp;         % 插值后采样率

rng(64);

%% ================= 多用户真实 DOA 设置 =================
% 注意：DOA 不是发送端调制的一部分，而是后续接收端阵列流形的参数。
% 这里统一定义真实角度，方便后续接收端仿真和估计误差对比。

user_angles_deg = [
     20,  10;     % User 1: theta, phi
    -15,   5;     % User 2: theta, phi
     35, -12      % User 3: theta, phi
];

theta_users_true = user_angles_deg(:,1) * pi/180;
phi_users_true   = user_angles_deg(:,2) * pi/180;

%% ================= 多用户导频发送端主流程 =================

% 1. 生成多用户已知导频比特
pilot_bits_users = generate_multiuser_pilot_bits( ...
    U, N_pilot, bits_per_sym, "random");

% 2. QPSK 导频映射
pilot_symbols_users = qpsk_modulate_multiuser_bits(pilot_bits_users);

% 3. 生成 RRC 成形滤波器
rrc = generate_rrc_filter(rolloff, span, sps);

% 4. 每个用户导频符号做 RRC 成形，得到复基带导频波形
pilot_bb_users = pulse_shape_multiuser_symbols( ...
    pilot_symbols_users, rrc, sps);

% 5. 插值到高采样率
pilot_bb_hi_users = interpolate_multiuser_baseband( ...
    pilot_bb_users, interp);

% 6. 构造时间轴
t = generate_time_axis(size(pilot_bb_hi_users, 1), fs);

% 7. 每个用户上变频到实射频导频信号
pilot_rf_users = upconvert_multiuser_to_rf( ...
    pilot_bb_hi_users, fc, fs);

% 8. 仅用于观察：多用户导频波形直接叠加
pilot_bb_hi_sum = sum(pilot_bb_hi_users, 2);
pilot_rf_sum    = sum(pilot_rf_users, 2);

%% ================= 输出基本信息 =================

fprintf('\n========== Multiuser DOA Pilot Transmitter ==========\n');
fprintf('Number of users              = %d\n', U);
fprintf('Pilot symbols per user       = %d\n', N_pilot);
fprintf('Bits per user                = %d\n', N_pilot * bits_per_sym);
fprintf('Symbol rate                  = %.2f MHz\n', Rs/1e6);
fprintf('Baseband sample rate         = %.2f MHz\n', fs_bb/1e6);
fprintf('RF simulation sample rate    = %.2f MHz\n', fs/1e6);
fprintf('Carrier frequency            = %.2f MHz\n', fc/1e6);
fprintf('Length of pilot_bb_hi_users  = %d samples\n', size(pilot_bb_hi_users, 1));

for u = 1:U
    fprintf('User %d true theta = %.2f deg, true phi = %.2f deg\n', ...
        u, user_angles_deg(u,1), user_angles_deg(u,2));
end

%% ================= 发送端可视化 =================

plot_multiuser_pilot_transmitter_waveforms( ...
    t, pilot_bb_hi_users, pilot_rf_users, ...
    pilot_bb_hi_sum, pilot_rf_sum, pilot_symbols_users, fc);



%% ================= 多用户 DOA 接收端参数设置 =================

RIS_row = 8;
RIS_col = 8;
M_ris = RIS_row * RIS_col;

lambda = 1;
d = 0.5 * lambda;

ref_amp = 1.5;

lpf_cutoff = 5e6;
fir_order = 800;

% 角度搜索范围
theta_grid_deg = -60:2:60;
phi_grid_deg   = -30:2:30;

theta_grid = theta_grid_deg * pi/180;
phi_grid   = phi_grid_deg * pi/180;

% 交替网格搜索迭代次数
alt_max_iter = 3;

% 局部连续优化最大迭代次数
local_max_iter = 100;

%% ================= 生成真实多用户阵列流形 =================

H_true = generate_multiuser_ris_steering_matrix( ...
    RIS_row, RIS_col, theta_users_true, phi_users_true, d, lambda);

%% ================= 生成 RIS 本地参考信号 =================

ref_symbols_ris = generate_reference_symbols_ris( ...
    N_pilot, M_ris, ref_amp);

%% ================= 计算符号能量窗口 =================

[sym_center_idx, win_len] = get_symbol_windows( ...
    N_pilot, sps, interp, span);

%% ================= 仿真得到接收端符号级能量观测 =================

E_obs = simulate_multiuser_pilot_symbol_energy_receiver( ...
    pilot_bb_hi_users, H_true, ref_symbols_ris, rrc, sps, interp, ...
    fc, fs, sym_center_idx, win_len, lpf_cutoff, fir_order);

%% ================= 多用户 DOA 估计算法对比 =================

% 预生成每个 RIS 单元的参考基带波形，三个算法复用，避免重复计算
ref_bb_hi_mat = generate_reference_waveforms_ris( ...
    ref_symbols_ris, rrc, sps, interp, size(pilot_bb_hi_users,1));

% ---------- 算法1：逐用户独立网格搜索 baseline ----------
[theta_hat_ind, phi_hat_ind, metric_ind] = estimate_multiuser_doa_independent_grid( ...
    E_obs, pilot_bb_hi_users, ref_bb_hi_mat, ...
    RIS_row, RIS_col, d, lambda, ...
    sym_center_idx, win_len, theta_grid, phi_grid);

% ---------- 算法2：多用户交替网格 ML ----------
[theta_hat_alt, phi_hat_alt, metric_hist_alt] = estimate_multiuser_doa_alternating_grid_ml( ...
    E_obs, pilot_bb_hi_users, ref_bb_hi_mat, ...
    RIS_row, RIS_col, d, lambda, ...
    sym_center_idx, win_len, theta_grid, phi_grid, ...
    theta_hat_ind, phi_hat_ind, alt_max_iter);

% ---------- 算法3：基于算法2结果的连续局部 ML 精修 ----------
[theta_hat_refine, phi_hat_refine, metric_refine] = refine_multiuser_doa_local_ml( ...
    E_obs, pilot_bb_hi_users, ref_bb_hi_mat, ...
    RIS_row, RIS_col, d, lambda, ...
    sym_center_idx, win_len, ...
    theta_hat_alt, phi_hat_alt, ...
    min(theta_grid), max(theta_grid), ...
    min(phi_grid), max(phi_grid), ...
    local_max_iter);

%% ================= GS / GN 两阶段 DOA 估计对比 =================
% 思路：
%   Step 1: 由符号级能量观测恢复多用户阵列响应矩阵 H_ris
%   Step 2: 将每个用户的恢复阵列响应拟合到 RIS 阵列流形，得到 DOA

Cp = estimate_pulse_energy_coefficient( ...
    N_pilot, rrc, sps, interp, sym_center_idx, win_len);

% ---------- GS-based channel recovery + DOA fitting ----------
t0_gs = 1000;

[theta_hat_gs, phi_hat_gs, H_hat_gs] = estimate_multiuser_doa_gs_channel_recovery( ...
    E_obs, pilot_symbols_users, ref_symbols_ris, Cp, ...
    RIS_row, RIS_col, d, lambda, theta_grid, phi_grid, t0_gs);

% ---------- GN-based channel recovery + DOA fitting ----------
maxIter_gn = 300;

[theta_hat_gn, phi_hat_gn, H_hat_gn] = estimate_multiuser_doa_gn_channel_recovery( ...
    E_obs, pilot_symbols_users, ref_symbols_ris, Cp, ...
    RIS_row, RIS_col, d, lambda, theta_grid, phi_grid, maxIter_gn);

%% 扩展打印结果
fprintf('\n========== GS / GN DOA Estimation Results ==========\n');

for u = 1:U
    fprintf('\nUser %d\n', u);
    fprintf('True : theta = %8.3f deg, phi = %8.3f deg\n', ...
        theta_users_true(u)*180/pi, phi_users_true(u)*180/pi);

    fprintf('GS   : theta = %8.3f deg, phi = %8.3f deg, err = (%+.3f, %+.3f) deg\n', ...
        theta_hat_gs(u)*180/pi, phi_hat_gs(u)*180/pi, ...
        (theta_hat_gs(u)-theta_users_true(u))*180/pi, ...
        (phi_hat_gs(u)-phi_users_true(u))*180/pi);

    fprintf('GN   : theta = %8.3f deg, phi = %8.3f deg, err = (%+.3f, %+.3f) deg\n', ...
        theta_hat_gn(u)*180/pi, phi_hat_gn(u)*180/pi, ...
        (theta_hat_gn(u)-theta_users_true(u))*180/pi, ...
        (phi_hat_gn(u)-phi_users_true(u))*180/pi);
end

fprintf('\nAverage GS abs error: theta %.3f deg, phi %.3f deg\n', ...
    mean(abs(theta_hat_gs-theta_users_true))*180/pi, ...
    mean(abs(phi_hat_gs-phi_users_true))*180/pi);

fprintf('Average GN abs error: theta %.3f deg, phi %.3f deg\n', ...
    mean(abs(theta_hat_gn-theta_users_true))*180/pi, ...
    mean(abs(phi_hat_gn-phi_users_true))*180/pi);

%% ================= 打印估计结果 =================

print_multiuser_doa_results( ...
    theta_users_true, phi_users_true, ...
    theta_hat_ind, phi_hat_ind, ...
    theta_hat_alt, phi_hat_alt, ...
    theta_hat_refine, phi_hat_refine);


% 增加打印
fprintf('\n========== Average DOA Error Summary ==========\n');

fprintf('Independent Grid : theta %.3f deg, phi %.3f deg\n', ...
    mean(abs(theta_hat_ind-theta_users_true))*180/pi, ...
    mean(abs(phi_hat_ind-phi_users_true))*180/pi);

fprintf('Alt. Grid ML     : theta %.3f deg, phi %.3f deg\n', ...
    mean(abs(theta_hat_alt-theta_users_true))*180/pi, ...
    mean(abs(phi_hat_alt-phi_users_true))*180/pi);

fprintf('Local ML Refine  : theta %.3f deg, phi %.3f deg\n', ...
    mean(abs(theta_hat_refine-theta_users_true))*180/pi, ...
    mean(abs(phi_hat_refine-phi_users_true))*180/pi);

fprintf('GS Channel Rec.  : theta %.3f deg, phi %.3f deg\n', ...
    mean(abs(theta_hat_gs-theta_users_true))*180/pi, ...
    mean(abs(phi_hat_gs-phi_users_true))*180/pi);

fprintf('GN Channel Rec.  : theta %.3f deg, phi %.3f deg\n', ...
    mean(abs(theta_hat_gn-theta_users_true))*180/pi, ...
    mean(abs(phi_hat_gn-phi_users_true))*180/pi);
%% ================= 可视化结果 =================

figure;
user_idx = 1:U;

subplot(2,1,1);
plot(user_idx, theta_users_true*180/pi, 'ko-', 'LineWidth', 1.5); hold on;
plot(user_idx, theta_hat_ind*180/pi, 'rx--', 'LineWidth', 1.2);
plot(user_idx, theta_hat_alt*180/pi, 'bs--', 'LineWidth', 1.2);
plot(user_idx, theta_hat_refine*180/pi, 'md--', 'LineWidth', 1.2);
grid on;
xlabel('User index');
ylabel('\theta / deg');
legend('True', 'Independent Grid', 'Alternating Grid ML', 'Local ML Refinement');
title('Multiuser DOA Estimation: \theta');

subplot(2,1,2);
plot(user_idx, phi_users_true*180/pi, 'ko-', 'LineWidth', 1.5); hold on;
plot(user_idx, phi_hat_ind*180/pi, 'rx--', 'LineWidth', 1.2);
plot(user_idx, phi_hat_alt*180/pi, 'bs--', 'LineWidth', 1.2);
plot(user_idx, phi_hat_refine*180/pi, 'md--', 'LineWidth', 1.2);
grid on;
xlabel('User index');
ylabel('\phi / deg');
legend('True', 'Independent Grid', 'Alternating Grid ML', 'Local ML Refinement');
title('Multiuser DOA Estimation: \phi');

figure;
plot(metric_hist_alt, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Alternating iteration');
ylabel('Normalized residual metric');
title('Alternating Grid ML Metric History');

figure;
user_idx = 1:U;

subplot(2,1,1);
plot(user_idx, theta_users_true*180/pi, 'ko-', 'LineWidth', 1.5); hold on;
plot(user_idx, theta_hat_alt*180/pi, 'bs--', 'LineWidth', 1.2);
plot(user_idx, theta_hat_refine*180/pi, 'md--', 'LineWidth', 1.2);
plot(user_idx, theta_hat_gs*180/pi, 'g^--', 'LineWidth', 1.2);
plot(user_idx, theta_hat_gn*180/pi, 'cv--', 'LineWidth', 1.2);
grid on;
xlabel('User index');
ylabel('\theta / deg');
legend('True', 'Alt. Grid ML', 'Local ML', 'GS-channel', 'GN-channel');
title('Multiuser DOA Estimation Comparison: \theta');

subplot(2,1,2);
plot(user_idx, phi_users_true*180/pi, 'ko-', 'LineWidth', 1.5); hold on;
plot(user_idx, phi_hat_alt*180/pi, 'bs--', 'LineWidth', 1.2);
plot(user_idx, phi_hat_refine*180/pi, 'md--', 'LineWidth', 1.2);
plot(user_idx, phi_hat_gs*180/pi, 'g^--', 'LineWidth', 1.2);
plot(user_idx, phi_hat_gn*180/pi, 'cv--', 'LineWidth', 1.2);
grid on;
xlabel('User index');
ylabel('\phi / deg');
legend('True', 'Alt. Grid ML', 'Local ML', 'GS-channel', 'GN-channel');
title('Multiuser DOA Estimation Comparison: \phi');


%% ========================================================================
%                    多用户 DOA 导频发送端函数区
% ========================================================================

function bits_users = generate_multiuser_pilot_bits(U, N_pilot, bits_per_sym, mode)
%GENERATE_MULTIUSER_PILOT_BITS 生成多用户导频比特
%
% 输入：
%   U            : 用户数
%   N_pilot      : 每个用户导频符号数
%   bits_per_sym : 每个符号对应比特数
%   mode         : "random" 或 "orthogonal_like"
%
% 输出：
%   bits_users   : N_bits x U，每一列对应一个用户的导频比特
%
% 说明：
%   - "random"：每个用户生成随机导频比特，但接收端已知；
%   - "orthogonal_like"：生成简单循环移位导频，用于增强多用户区分度。
%     注意这里不是严格正交 QPSK 序列，只是提供一个结构化导频示例。

    N_bits = N_pilot * bits_per_sym;

    if nargin < 4
        mode = "random";
    end

    bits_users = zeros(N_bits, U);

    switch mode

        case "random"
            bits_users = randi([0 1], N_bits, U);

        case "orthogonal_like"

            base_bits = randi([0 1], N_bits, 1);

            for u = 1:U
                shift_len = round((u-1) * N_bits / U);
                bits_users(:,u) = circshift(base_bits, shift_len);
            end

        otherwise
            error('未知导频比特生成模式：%s', mode);
    end
end


function symbols_users = qpsk_modulate_multiuser_bits(bits_users)
%QPSK_MODULATE_MULTIUSER_BITS 多用户 QPSK 映射
%
% 输入：
%   bits_users    : N_bits x U，每一列对应一个用户
%
% 输出：
%   symbols_users : N_sym x U，每一列对应一个用户的 QPSK 符号

    [N_bits, U] = size(bits_users);

    if mod(N_bits, 2) ~= 0
        error('QPSK 调制要求每个用户的比特数为 2 的整数倍。');
    end

    N_sym = N_bits / 2;

    symbols_users = zeros(N_sym, U);

    for u = 1:U

        bits_u = bits_users(:,u);

        bits_reshape = reshape(bits_u, 2, []).';

        symbol_idx = bi2de(bits_reshape, 'left-msb');

        symbols_u = qammod(symbol_idx, 4, 'gray', ...
            'UnitAveragePower', true);

        symbols_users(:,u) = symbols_u(:);
    end
end


function rrc = generate_rrc_filter(rolloff, span, sps)
%GENERATE_RRC_FILTER 生成根升余弦成形滤波器
%
% 输入：
%   rolloff : 滚降系数
%   span    : 滤波器跨度，单位：符号
%   sps     : 每符号采样点数
%
% 输出：
%   rrc     : RRC 滤波器系数

    rrc = rcosdesign(rolloff, span, sps, 'sqrt');

    rrc = rrc(:);
end


function tx_bb_users = pulse_shape_multiuser_symbols(symbols_users, rrc, sps)
%PULSE_SHAPE_MULTIUSER_SYMBOLS 对每个用户符号做 RRC 成形
%
% 输入：
%   symbols_users : N_sym x U，多用户符号矩阵
%   rrc           : RRC 滤波器
%   sps           : 每符号采样点数
%
% 输出：
%   tx_bb_users   : N_sample x U，每列为一个用户的复基带波形

    [~, U] = size(symbols_users);

    tx_bb_1 = upfirdn(symbols_users(:,1), rrc, sps, 1);
    tx_bb_1 = tx_bb_1(:);

    N_sample = length(tx_bb_1);

    tx_bb_users = zeros(N_sample, U);

    tx_bb_users(:,1) = tx_bb_1;

    for u = 2:U
        tx_bb_u = upfirdn(symbols_users(:,u), rrc, sps, 1);
        tx_bb_users(:,u) = tx_bb_u(:);
    end
end


function tx_bb_hi_users = interpolate_multiuser_baseband(tx_bb_users, interp)
%INTERPOLATE_MULTIUSER_BASEBAND 多用户基带波形插值
%
% 输入：
%   tx_bb_users    : N_sample x U，多用户复基带波形
%   interp         : 插值倍数
%
% 输出：
%   tx_bb_hi_users : N_sample_hi x U，插值后的多用户复基带波形

    [~, U] = size(tx_bb_users);

    if interp == 1
        tx_bb_hi_users = tx_bb_users;
        return;
    end

    tx_bb_hi_1 = resample(tx_bb_users(:,1), interp, 1);
    tx_bb_hi_1 = tx_bb_hi_1(:);

    N_sample_hi = length(tx_bb_hi_1);

    tx_bb_hi_users = zeros(N_sample_hi, U);

    tx_bb_hi_users(:,1) = tx_bb_hi_1;

    for u = 2:U
        tx_bb_hi_u = resample(tx_bb_users(:,u), interp, 1);
        tx_bb_hi_users(:,u) = tx_bb_hi_u(:);
    end
end


function t = generate_time_axis(N_sample, fs)
%GENERATE_TIME_AXIS 生成时间轴
%
% 输入：
%   N_sample : 采样点数
%   fs       : 采样率
%
% 输出：
%   t        : N_sample x 1 时间向量

    t = (0:N_sample-1).' / fs;
end


function tx_rf_users = upconvert_multiuser_to_rf(tx_bb_hi_users, fc, fs)
%UPCONVERT_MULTIUSER_TO_RF 多用户复基带上变频为实射频信号
%
% 对每个用户：
%   x_RF,u(t) = Re{ x_BB,u(t) exp(j2πf_ct) }
%
% 输入：
%   tx_bb_hi_users : N_sample x U，多用户高采样率复基带波形
%   fc             : 载波频率
%   fs             : 采样率
%
% 输出：
%   tx_rf_users    : N_sample x U，多用户实射频信号

    [N_sample, U] = size(tx_bb_hi_users);

    t = (0:N_sample-1).' / fs;

    carrier = exp(1j * 2*pi*fc*t);

    tx_rf_users = zeros(N_sample, U);

    for u = 1:U
        tx_rf_users(:,u) = real(tx_bb_hi_users(:,u) .* carrier);
    end
end


function plot_multiuser_pilot_transmitter_waveforms( ...
    t, pilot_bb_hi_users, pilot_rf_users, ...
    pilot_bb_hi_sum, pilot_rf_sum, pilot_symbols_users, fc)
%PLOT_MULTIUSER_PILOT_TRANSMITTER_WAVEFORMS
% 绘制多用户 DOA 导频发送端关键波形

    [~, U] = size(pilot_bb_hi_users);

    figure('Position', [100, 100, 1200, 900]);

    subplot(5,1,1);
    hold on;
    for u = 1:U
        plot(real(pilot_symbols_users(:,u)), imag(pilot_symbols_users(:,u)), ...
            'o', 'LineWidth', 1.2);
    end
    grid on; axis equal;
    xlabel('In-phase');
    ylabel('Quadrature');
    title('Multiuser QPSK Pilot Symbols');
    legend(compose('User %d', 1:U));

    subplot(5,1,2);
    hold on;
    for u = 1:U
        plot(t*1e6, real(pilot_bb_hi_users(:,u)), 'LineWidth', 1.0);
    end
    grid on;
    xlabel('Time / \mus');
    ylabel('Amplitude');
    title('Real Part of User Pilot Baseband Signals');
    legend(compose('User %d', 1:U));

    subplot(5,1,3);
    hold on;
    for u = 1:U
        plot(t*1e6, abs(pilot_bb_hi_users(:,u)), 'LineWidth', 1.0);
    end
    grid on;
    xlabel('Time / \mus');
    ylabel('Envelope');
    title('Envelope of User Pilot Baseband Signals');
    legend(compose('User %d', 1:U));

    subplot(5,1,4);
    plot(t*1e6, real(pilot_bb_hi_sum), 'LineWidth', 1.1); hold on;
    plot(t*1e6, imag(pilot_bb_hi_sum), 'LineWidth', 1.1);
    grid on;
    xlabel('Time / \mus');
    ylabel('Amplitude');
    legend('I', 'Q');
    title('Direct Sum of Multiuser Pilot Complex Baseband Signals');

    subplot(5,1,5);

    if length(t) > 1
        T_show = max(2e-6, t(end));
        N_show = find(t <= T_show, 1, 'last');
        if isempty(N_show)
            N_show = min(length(t), 2000);
        end
    else
        N_show = 1;
    end

    plot(t(1:N_show)*1e6, pilot_rf_sum(1:N_show), 'LineWidth', 1.0);
    grid on;
    xlabel('Time / \mus');
    ylabel('Amplitude');
    title('Direct Sum of Multiuser Pilot RF Signals');

    sgtitle('Multiuser DOA Pilot Transmitter');
end


%% ========================================================================
%                    多用户 DOA 接收端函数区
% ========================================================================

function H_ris = generate_multiuser_ris_steering_matrix( ...
    RIS_row, RIS_col, theta_users, phi_users, d, lambda)
% 生成多用户 RIS 阵列流形矩阵
%
% 输出：
%   H_ris : M_ris x U

    theta_users = theta_users(:);
    phi_users = phi_users(:);

    U = length(theta_users);

    if length(phi_users) ~= U
        error('theta_users 和 phi_users 长度必须一致。');
    end

    M_ris = RIS_row * RIS_col;
    H_ris = zeros(M_ris, U);

    for u = 1:U
        H_ris(:,u) = generate_ris_steering_vector( ...
            RIS_row, RIS_col, theta_users(u), phi_users(u), d, lambda);
    end
end


function h_ris = generate_ris_steering_vector(RIS_row, RIS_col, theta, phi, d, lambda)
% 生成二维 RIS 平面阵列流形
%
% 假设 RIS 位于 y-z 平面，x 轴为法向。
% row 对应 z 方向，col 对应 y 方向。

    k0 = 2*pi/lambda;

    [y_idx, z_idx] = meshgrid(0:RIS_col-1, 0:RIS_row-1);

    y_pos = (y_idx(:) - (RIS_col-1)/2) * d;
    z_pos = (z_idx(:) - (RIS_row-1)/2) * d;

    ky = k0 * sin(theta) * cos(phi);
    kz = k0 * sin(theta) * sin(phi);

    h_ris = exp(1j * (ky * y_pos + kz * z_pos));

    h_ris = h_ris(:);
end


function ref_symbols_ris = generate_reference_symbols_ris(N_sym, M_ris, ref_amp)
% 为每个 RIS 单元生成已知参考符号序列

    k = (0:N_sym-1).';

    base_phase = mod(2*pi*0.137*k.^2 + pi/7*k, 2*pi);

    rng(2028);
    ris_phase = 2*pi*rand(1, M_ris);

    ref_symbols_ris = zeros(N_sym, M_ris);

    for m = 1:M_ris
        ref_symbols_ris(:,m) = ref_amp * exp(1j * (base_phase + ris_phase(m)));
    end
end


function [sym_center_idx, win_len] = get_symbol_windows(N_sym, sps, interp, span)
% 计算每个符号的中心采样点和积分窗口长度

    sym_samp_hi = sps * interp;

    rrc_delay_bb = span * sps / 2;
    rrc_delay_hi = rrc_delay_bb * interp;

    sym_center_idx = round(rrc_delay_hi + 1 + (0:N_sym-1).' * sym_samp_hi);

    win_len = sym_samp_hi;
end


function E_simo = simulate_multiuser_pilot_symbol_energy_receiver( ...
    pilot_bb_hi_users, H_ris, ref_symbols_ris, rrc, sps, interp, ...
    fc, fs, sym_center_idx, win_len, lpf_cutoff, fir_order)
% 仿真多用户导频阶段的 SIMO 符号级能量接收
%
% 对每个 RIS 单元：
%   x_m(t) = sum_u h_{m,u} x_u(t)
%   y_m(t) = x_m(t) + r_m(t)
%   E_{m,k} = symbol average of |y_m(t)|^2

    [N_sample, U] = size(pilot_bb_hi_users);
    [M_ris, U_h] = size(H_ris);
    [N_sym, M_ref] = size(ref_symbols_ris);

    if U_h ~= U
        error('H_ris 的列数必须等于用户数。');
    end

    if M_ref ~= M_ris
        error('ref_symbols_ris 的列数必须等于 RIS 单元数。');
    end

    E_simo = zeros(M_ris, N_sym);

    t = (0:N_sample-1).' / fs;
    carrier = exp(1j * 2*pi*fc*t);

    for m = 1:M_ris

        rx_unknown_bb_m = zeros(N_sample, 1);

        for u = 1:U
            rx_unknown_bb_m = rx_unknown_bb_m ...
                + H_ris(m,u) * pilot_bb_hi_users(:,u);
        end

        ref_bb_m = pulse_shape_symbols(ref_symbols_ris(:,m), rrc, sps);
        ref_bb_hi_m = interpolate_baseband(ref_bb_m, interp);

        L = min([length(rx_unknown_bb_m), length(ref_bb_hi_m), length(carrier)]);

        rx_unknown_bb_m = rx_unknown_bb_m(1:L);
        ref_bb_hi_m = ref_bb_hi_m(1:L);
        carrier_m = carrier(1:L);

        rx_unknown_rf_m = real(rx_unknown_bb_m .* carrier_m);
        ref_rf_m = real(ref_bb_hi_m .* carrier_m);

        mixed_rf_m = rx_unknown_rf_m + ref_rf_m;

        energy_waveform_m = square_law_detector( ...
            mixed_rf_m, fs, lpf_cutoff, fir_order);

        E_simo(m,:) = extract_symbol_energy_from_waveform( ...
            energy_waveform_m, sym_center_idx, win_len);
    end
end


function ref_bb_hi_mat = generate_reference_waveforms_ris( ...
    ref_symbols_ris, rrc, sps, interp, target_len)
% 生成每个 RIS 单元的参考复基带波形

    [~, M_ris] = size(ref_symbols_ris);

    ref_bb_hi_mat = zeros(target_len, M_ris);

    for m = 1:M_ris

        ref_bb_m = pulse_shape_symbols(ref_symbols_ris(:,m), rrc, sps);
        ref_bb_hi_m = interpolate_baseband(ref_bb_m, interp);

        L = min(target_len, length(ref_bb_hi_m));

        ref_bb_hi_mat(1:L,m) = ref_bb_hi_m(1:L);
    end
end


function [theta_hat, phi_hat, metric_maps] = estimate_multiuser_doa_independent_grid( ...
    E_obs, pilot_bb_hi_users, ref_bb_hi_mat, ...
    RIS_row, RIS_col, d, lambda, ...
    sym_center_idx, win_len, theta_grid, phi_grid)
% 算法1：逐用户独立二维网格搜索 baseline
%
% 对每个用户 u，忽略其它用户，搜索使得：
%   E_obs ≈ E_pred(h_u, x_u, ref)
% 残差最小的角度。
%
% 这个算法不是严格多用户 ML，只作为 baseline。

    [~, U] = size(pilot_bb_hi_users);
    [M_ris, ~] = size(E_obs);

    theta_hat = zeros(U, 1);
    phi_hat = zeros(U, 1);

    metric_maps = cell(U, 1);

    for u = 1:U

        metric_map_u = zeros(length(theta_grid), length(phi_grid));

        pilot_single = pilot_bb_hi_users(:,u);

        for it = 1:length(theta_grid)

            theta = theta_grid(it);

            for ip = 1:length(phi_grid)

                phi = phi_grid(ip);

                h_cand = generate_ris_steering_vector( ...
                    RIS_row, RIS_col, theta, phi, d, lambda);

                H_cand = h_cand;   % M_ris x 1

                E_pred = predict_multiuser_symbol_energy_baseband( ...
                    pilot_single, H_cand, ref_bb_hi_mat, sym_center_idx, win_len);

                residual = E_obs - E_pred;

                metric_map_u(it, ip) = norm(residual, 'fro')^2 ...
                    / max(norm(E_obs, 'fro')^2, 1e-12);
            end
        end

        [~, idx_min] = min(metric_map_u(:));
        [idx_theta, idx_phi] = ind2sub(size(metric_map_u), idx_min);

        theta_hat(u) = theta_grid(idx_theta);
        phi_hat(u) = phi_grid(idx_phi);

        metric_maps{u} = metric_map_u;
    end
end


function [theta_hat, phi_hat, metric_hist] = estimate_multiuser_doa_alternating_grid_ml( ...
    E_obs, pilot_bb_hi_users, ref_bb_hi_mat, ...
    RIS_row, RIS_col, d, lambda, ...
    sym_center_idx, win_len, theta_grid, phi_grid, ...
    theta_init, phi_init, max_iter)
% 算法2：多用户交替网格 ML
%
% 初始化为独立搜索结果。
% 每次固定其它用户角度，只对当前用户做二维网格搜索。
% 目标函数是完整多用户能量残差。

    [~, U] = size(pilot_bb_hi_users);

    theta_hat = theta_init(:);
    phi_hat = phi_init(:);

    metric_hist = zeros(max_iter, 1);

    for iter = 1:max_iter

        for u = 1:U

            best_metric = inf;
            best_theta = theta_hat(u);
            best_phi = phi_hat(u);

            theta_temp = theta_hat;
            phi_temp = phi_hat;

            for it = 1:length(theta_grid)

                theta_temp(u) = theta_grid(it);

                for ip = 1:length(phi_grid)

                    phi_temp(u) = phi_grid(ip);

                    metric_val = multiuser_doa_metric( ...
                        E_obs, pilot_bb_hi_users, ref_bb_hi_mat, ...
                        RIS_row, RIS_col, d, lambda, ...
                        sym_center_idx, win_len, theta_temp, phi_temp);

                    if metric_val < best_metric
                        best_metric = metric_val;
                        best_theta = theta_temp(u);
                        best_phi = phi_temp(u);
                    end
                end
            end

            theta_hat(u) = best_theta;
            phi_hat(u) = best_phi;
        end

        metric_hist(iter) = multiuser_doa_metric( ...
            E_obs, pilot_bb_hi_users, ref_bb_hi_mat, ...
            RIS_row, RIS_col, d, lambda, ...
            sym_center_idx, win_len, theta_hat, phi_hat);

        fprintf('Alternating ML iter %d/%d, metric = %.4e\n', ...
            iter, max_iter, metric_hist(iter));
    end
end


function [theta_hat, phi_hat, metric_final] = refine_multiuser_doa_local_ml( ...
    E_obs, pilot_bb_hi_users, ref_bb_hi_mat, ...
    RIS_row, RIS_col, d, lambda, ...
    sym_center_idx, win_len, ...
    theta_init, phi_init, ...
    theta_min, theta_max, phi_min, phi_max, max_iter)
% 算法3：连续局部 ML 精修
%
% 使用 fminsearch 从交替网格结果出发，优化所有用户角度。
% 这里加入越界惩罚，使搜索限制在给定角度范围内。

    theta_init = theta_init(:);
    phi_init = phi_init(:);

    U = length(theta_init);

    x0 = zeros(2*U, 1);
    x0(1:U) = theta_init;
    x0(U+1:end) = phi_init;

    obj_fun = @(x) local_doa_objective_with_bounds( ...
        x, E_obs, pilot_bb_hi_users, ref_bb_hi_mat, ...
        RIS_row, RIS_col, d, lambda, ...
        sym_center_idx, win_len, ...
        theta_min, theta_max, phi_min, phi_max);

    options = optimset( ...
        'MaxIter', max_iter, ...
        'MaxFunEvals', 5e4, ...
        'Display', 'iter');

    x_hat = fminsearch(obj_fun, x0, options);

    theta_hat = x_hat(1:U);
    phi_hat = x_hat(U+1:end);

    theta_hat = min(max(theta_hat, theta_min), theta_max);
    phi_hat = min(max(phi_hat, phi_min), phi_max);

    metric_final = multiuser_doa_metric( ...
        E_obs, pilot_bb_hi_users, ref_bb_hi_mat, ...
        RIS_row, RIS_col, d, lambda, ...
        sym_center_idx, win_len, theta_hat, phi_hat);
end


function obj = local_doa_objective_with_bounds( ...
    x, E_obs, pilot_bb_hi_users, ref_bb_hi_mat, ...
    RIS_row, RIS_col, d, lambda, ...
    sym_center_idx, win_len, ...
    theta_min, theta_max, phi_min, phi_max)
% 带边界惩罚的连续 DOA 目标函数

    U = length(x) / 2;

    theta = x(1:U);
    phi = x(U+1:end);

    penalty = 0;

    if any(theta < theta_min)
        penalty = penalty + 1e3 * sum((theta(theta < theta_min) - theta_min).^2);
    end

    if any(theta > theta_max)
        penalty = penalty + 1e3 * sum((theta(theta > theta_max) - theta_max).^2);
    end

    if any(phi < phi_min)
        penalty = penalty + 1e3 * sum((phi(phi < phi_min) - phi_min).^2);
    end

    if any(phi > phi_max)
        penalty = penalty + 1e3 * sum((phi(phi > phi_max) - phi_max).^2);
    end

    theta = min(max(theta, theta_min), theta_max);
    phi = min(max(phi, phi_min), phi_max);

    metric_val = multiuser_doa_metric( ...
        E_obs, pilot_bb_hi_users, ref_bb_hi_mat, ...
        RIS_row, RIS_col, d, lambda, ...
        sym_center_idx, win_len, theta, phi);

    obj = metric_val + penalty;
end


function metric_val = multiuser_doa_metric( ...
    E_obs, pilot_bb_hi_users, ref_bb_hi_mat, ...
    RIS_row, RIS_col, d, lambda, ...
    sym_center_idx, win_len, theta_users, phi_users)
% 多用户 DOA 估计残差目标函数

    H_cand = generate_multiuser_ris_steering_matrix( ...
        RIS_row, RIS_col, theta_users, phi_users, d, lambda);

    E_pred = predict_multiuser_symbol_energy_baseband( ...
        pilot_bb_hi_users, H_cand, ref_bb_hi_mat, sym_center_idx, win_len);

    residual = E_obs - E_pred;

    metric_val = norm(residual, 'fro')^2 / max(norm(E_obs, 'fro')^2, 1e-12);
end


function E_pred = predict_multiuser_symbol_energy_baseband( ...
    pilot_bb_hi_users, H_ris, ref_bb_hi_mat, sym_center_idx, win_len)
% 根据候选多用户阵列流形预测符号级能量
%
% 模型：
%   u_m[n] = sum_u h_{m,u} x_u[n] + r_m[n]
%   E_{m,k} = mean_{n in window k} |u_m[n]|^2

    if isvector(pilot_bb_hi_users)
        pilot_bb_hi_users = pilot_bb_hi_users(:);
    end

    [N_sample, U] = size(pilot_bb_hi_users);
    [M_ris, U_h] = size(H_ris);

    if U_h ~= U
        error('H_ris 的列数必须等于 pilot_bb_hi_users 的用户数。');
    end

    L = min(N_sample, size(ref_bb_hi_mat,1));

    pilot_bb_hi_users = pilot_bb_hi_users(1:L,:);
    ref_bb_hi_mat = ref_bb_hi_mat(1:L,:);

    N_sym = length(sym_center_idx);

    E_pred = zeros(M_ris, N_sym);

    for m = 1:M_ris

        u_m = zeros(L, 1);

        for u = 1:U
            u_m = u_m + H_ris(m,u) * pilot_bb_hi_users(:,u);
        end

        u_m = u_m + ref_bb_hi_mat(:,m);

        energy_m = abs(u_m).^2;

        E_pred(m,:) = extract_symbol_energy_from_waveform( ...
            energy_m, sym_center_idx, win_len);
    end
end


function energy_est = square_law_detector(rx_rf, fs, lpf_cutoff, fir_order)
% 平方律检波器

    rx_rf = rx_rf(:);

    rx_square = rx_rf.^2;

    b_lpf = fir1(fir_order, lpf_cutoff/(fs/2), 'low');

    rx_energy_lpf = filtfilt(b_lpf, 1, rx_square);

    energy_est = 2 * rx_energy_lpf;

    energy_est = max(real(energy_est), 0);

    energy_est = energy_est(:);
end


function E_symbol = extract_symbol_energy_from_waveform(energy_waveform, sym_center_idx, win_len)
% 从能量波形中提取每个符号平均能量

    energy_waveform = energy_waveform(:);

    N_sym = length(sym_center_idx);

    E_symbol = zeros(1, N_sym);

    for k = 1:N_sym

        cidx = sym_center_idx(k);

        left_idx  = round(cidx - win_len/2);
        right_idx = round(cidx + win_len/2 - 1);

        left_idx  = max(left_idx, 1);
        right_idx = min(right_idx, length(energy_waveform));

        if left_idx <= right_idx
            E_symbol(k) = mean(energy_waveform(left_idx:right_idx));
        else
            E_symbol(k) = NaN;
        end
    end
end


function tx_bb = pulse_shape_symbols(symbols, rrc, sps)
% 符号上采样并 RRC 成形

    symbols = symbols(:);

    tx_bb = upfirdn(symbols, rrc, sps, 1);

    tx_bb = tx_bb(:);
end


function tx_bb_hi = interpolate_baseband(tx_bb, interp)
% 基带信号插值

    tx_bb = tx_bb(:);

    if interp == 1
        tx_bb_hi = tx_bb;
    else
        tx_bb_hi = resample(tx_bb, interp, 1);
    end

    tx_bb_hi = tx_bb_hi(:);
end


function print_multiuser_doa_results( ...
    theta_true, phi_true, ...
    theta_ind, phi_ind, ...
    theta_alt, phi_alt, ...
    theta_refine, phi_refine)
% 打印多用户 DOA 估计结果

    theta_true = theta_true(:);
    phi_true = phi_true(:);

    U = length(theta_true);

    fprintf('\n========== Multiuser DOA Estimation Results ==========\n');

    for u = 1:U

        fprintf('\nUser %d\n', u);
        fprintf('True        : theta = %8.3f deg, phi = %8.3f deg\n', ...
            theta_true(u)*180/pi, phi_true(u)*180/pi);

        fprintf('Independent : theta = %8.3f deg, phi = %8.3f deg, ', ...
            theta_ind(u)*180/pi, phi_ind(u)*180/pi);
        fprintf('err = (%+.3f, %+.3f) deg\n', ...
            (theta_ind(u)-theta_true(u))*180/pi, ...
            (phi_ind(u)-phi_true(u))*180/pi);

        fprintf('Alt. Grid   : theta = %8.3f deg, phi = %8.3f deg, ', ...
            theta_alt(u)*180/pi, phi_alt(u)*180/pi);
        fprintf('err = (%+.3f, %+.3f) deg\n', ...
            (theta_alt(u)-theta_true(u))*180/pi, ...
            (phi_alt(u)-phi_true(u))*180/pi);

        fprintf('Refined ML  : theta = %8.3f deg, phi = %8.3f deg, ', ...
            theta_refine(u)*180/pi, phi_refine(u)*180/pi);
        fprintf('err = (%+.3f, %+.3f) deg\n', ...
            (theta_refine(u)-theta_true(u))*180/pi, ...
            (phi_refine(u)-phi_true(u))*180/pi);
    end

    fprintf('\nAverage absolute angle errors:\n');

    fprintf('Independent : theta %.3f deg, phi %.3f deg\n', ...
        mean(abs(theta_ind-theta_true))*180/pi, ...
        mean(abs(phi_ind-phi_true))*180/pi);

    fprintf('Alt. Grid   : theta %.3f deg, phi %.3f deg\n', ...
        mean(abs(theta_alt-theta_true))*180/pi, ...
        mean(abs(phi_alt-phi_true))*180/pi);

    fprintf('Refined ML  : theta %.3f deg, phi %.3f deg\n', ...
        mean(abs(theta_refine-theta_true))*180/pi, ...
        mean(abs(phi_refine-phi_true))*180/pi);
end


function [theta_hat, phi_hat, H_hat] = estimate_multiuser_doa_gs_channel_recovery( ...
    E_obs, pilot_symbols_users, ref_symbols_ris, Cp, ...
    RIS_row, RIS_col, d, lambda, theta_grid, phi_grid, t0)
%ESTIMATE_MULTIUSER_DOA_GS_CHANNEL_RECOVERY
% 使用 GS 先恢复多用户阵列响应矩阵 H_hat，再进行 DOA 网格拟合。
%
% 符号级幅值模型：
%   sqrt(E_{m,k}) ≈ | sqrt(Cp)*(sum_u h_{m,u}s_{u,k} + b_{m,k}) |
%
% 将未知 H_ris 展开成向量 g：
%   g = [h_1; h_2; ...; h_U]
%
% 构造：
%   z = |A*g + b|
%
% 如果旧 biased_gs_algorithm 使用：
%   z = abs(A_gs' * g + b)
%
% 则传入：
%   A_gs = A.'

    [A_pr, b_pr, z_pr, M_ris, U] = build_multiuser_channel_pr_model( ...
        E_obs, pilot_symbols_users, ref_symbols_ris, Cp);

    % 旧 GS 接口：z = abs(A_gs' * g + b)
    A_gs = A_pr.';

    g_hat = biased_gs_algorithm(z_pr, A_gs, b_pr, t0);

    g_hat = g_hat(:);

    H_hat = reshape_channel_vector_to_matrix(g_hat, M_ris, U);

    [theta_hat, phi_hat] = fit_multiuser_channel_to_doa_grid( ...
        H_hat, RIS_row, RIS_col, d, lambda, theta_grid, phi_grid);
end


function [theta_hat, phi_hat, H_hat] = estimate_multiuser_doa_gn_channel_recovery( ...
    E_obs, pilot_symbols_users, ref_symbols_ris, Cp, ...
    RIS_row, RIS_col, d, lambda, theta_grid, phi_grid, maxIter)
%ESTIMATE_MULTIUSER_DOA_GN_CHANNEL_RECOVERY
% 使用 GN 先恢复多用户阵列响应矩阵 H_hat，再进行 DOA 网格拟合。
%
% biased_gn_algorithm 的接口应为：
%   z = abs(A*g + b)

    [A_pr, b_pr, z_pr, M_ris, U] = build_multiuser_channel_pr_model( ...
        E_obs, pilot_symbols_users, ref_symbols_ris, Cp);

    g_hat = biased_gn_algorithm(z_pr, A_pr, b_pr, maxIter);

    g_hat = g_hat(:);

    H_hat = reshape_channel_vector_to_matrix(g_hat, M_ris, U);

    [theta_hat, phi_hat] = fit_multiuser_channel_to_doa_grid( ...
        H_hat, RIS_row, RIS_col, d, lambda, theta_grid, phi_grid);
end


function [A_pr, b_pr, z_pr, M_ris, U] = build_multiuser_channel_pr_model( ...
    E_obs, pilot_symbols_users, ref_symbols_ris, Cp)
%BUILD_MULTIUSER_CHANNEL_PR_MODEL
% 构造多用户 DOA 中用于 GS/GN 的阵列响应恢复模型。
%
% 观测模型：
%   E_{m,k} ≈ Cp * | sum_u h_{m,u}s_{u,k} + b_{m,k} |^2
%
% 转成幅值模型：
%   sqrt(E_{m,k}) ≈ | sqrt(Cp)*sum_u h_{m,u}s_{u,k}
%                    + sqrt(Cp)*b_{m,k} |
%
% 令：
%   g = [h_1; h_2; ...; h_U] ∈ C^{M_ris*U}
%
% 对每个观测 (m,k)：
%   z_{m,k} = | A(row,:) * g + b(row) |
%
% 其中 A(row, (u-1)*M_ris + m) = sqrt(Cp)*s_{k,u}

    [M_ris, N_sym] = size(E_obs);
    [N_sym_pilot, U] = size(pilot_symbols_users);

    if N_sym_pilot ~= N_sym
        error('pilot_symbols_users 的行数必须等于 E_obs 的列数。');
    end

    if size(ref_symbols_ris,1) ~= N_sym || size(ref_symbols_ris,2) ~= M_ris
        error('ref_symbols_ris 的尺寸必须为 N_sym x M_ris。');
    end

    if Cp <= 0
        error('Cp 必须为正数。');
    end

    N_obs = M_ris * N_sym;
    N_unknown = M_ris * U;

    A_pr = zeros(N_obs, N_unknown);
    b_pr = zeros(N_obs, 1);
    z_pr = zeros(N_obs, 1);

    row = 0;

    for k = 1:N_sym
        for m = 1:M_ris

            row = row + 1;

            % 幅值观测
            z_pr(row) = sqrt(max(E_obs(m,k), 0));

            % 已知参考偏置
            b_pr(row) = sqrt(Cp) * ref_symbols_ris(k,m);

            % 多用户阵列响应未知量的线性系数
            for u = 1:U
                col = (u-1)*M_ris + m;
                A_pr(row, col) = sqrt(Cp) * pilot_symbols_users(k,u);
            end
        end
    end
end

function H_hat = reshape_channel_vector_to_matrix(g_hat, M_ris, U)
%RESHAPE_CHANNEL_VECTOR_TO_MATRIX
% 将 g = [h_1; h_2; ...; h_U] 还原成 H_hat = [h_1, h_2, ..., h_U]

    g_hat = g_hat(:);

    if length(g_hat) ~= M_ris * U
        error('g_hat 长度必须等于 M_ris * U。');
    end

    H_hat = zeros(M_ris, U);

    for u = 1:U
        idx = (u-1)*M_ris + (1:M_ris);
        H_hat(:,u) = g_hat(idx);
    end
end


function [theta_hat, phi_hat] = fit_multiuser_channel_to_doa_grid( ...
    H_hat, RIS_row, RIS_col, d, lambda, theta_grid, phi_grid)
%FIT_MULTIUSER_CHANNEL_TO_DOA_GRID
% 将每个用户恢复出来的阵列响应 h_hat 拟合到阵列流形 a(theta,phi)。
%
% 对每个用户：
%   min_{theta,phi,alpha} || h_hat - alpha*a(theta,phi) ||^2
%
% 其中 alpha 是复数尺度因子，用于吸收幅度和全局相位误差。

    [M_ris, U] = size(H_hat);

    if M_ris ~= RIS_row * RIS_col
        error('H_hat 行数必须等于 RIS_row * RIS_col。');
    end

    theta_hat = zeros(U, 1);
    phi_hat = zeros(U, 1);

    for u = 1:U

        h_u = H_hat(:,u);

        [theta_hat(u), phi_hat(u)] = fit_single_channel_to_doa_grid( ...
            h_u, RIS_row, RIS_col, d, lambda, theta_grid, phi_grid);
    end
end

function [theta_hat, phi_hat, metric_map] = fit_single_channel_to_doa_grid( ...
    h_hat, RIS_row, RIS_col, d, lambda, theta_grid, phi_grid)
%FIT_SINGLE_CHANNEL_TO_DOA_GRID
% 将恢复出的阵列响应 h_hat 拟合到二维阵列流形。
%
% 目标：
%   min_{theta,phi,alpha} || h_hat - alpha*a(theta,phi) ||^2
%
% 对固定 a(theta,phi)，最优：
%   alpha = (a^H h_hat) / (a^H a)

    h_hat = h_hat(:);

    metric_map = zeros(length(theta_grid), length(phi_grid));

    best_metric = inf;
    theta_hat = theta_grid(1);
    phi_hat = phi_grid(1);

    norm_h = norm(h_hat)^2;

    if norm_h < 1e-12
        warning('h_hat 接近零向量，DOA 拟合可能不可靠。');
        return;
    end

    for it = 1:length(theta_grid)

        theta = theta_grid(it);

        for ip = 1:length(phi_grid)

            phi = phi_grid(ip);

            a = generate_ris_steering_vector( ...
                RIS_row, RIS_col, theta, phi, d, lambda);

            % 最优复数尺度因子
            alpha = (a' * h_hat) / max(a' * a, 1e-12);

            residual = h_hat - alpha * a;

            metric = norm(residual)^2 / max(norm_h, 1e-12);

            metric_map(it, ip) = metric;

            if metric < best_metric
                best_metric = metric;
                theta_hat = theta;
                phi_hat = phi;
            end
        end
    end
end

function Cp = estimate_pulse_energy_coefficient( ...
    N_sym, rrc, sps, interp, sym_center_idx, win_len)
%ESTIMATE_PULSE_ENERGY_COEFFICIENT
% 估计符号窗口内的脉冲平均能量系数。
%
% 用一个单位符号通过同样的 RRC 成形和插值，计算其在符号窗口内的平均能量。

    mid_k = ceil(N_sym/2);

    s_basis = zeros(N_sym, 1);
    s_basis(mid_k) = 1;

    bb_basis = pulse_shape_symbols(s_basis, rrc, sps);
    bb_basis_hi = interpolate_baseband(bb_basis, interp);

    cidx = sym_center_idx(mid_k);

    left_idx  = round(cidx - win_len/2);
    right_idx = round(cidx + win_len/2 - 1);

    left_idx  = max(left_idx, 1);
    right_idx = min(right_idx, length(bb_basis_hi));

    Cp = mean(abs(bb_basis_hi(left_idx:right_idx)).^2);
end