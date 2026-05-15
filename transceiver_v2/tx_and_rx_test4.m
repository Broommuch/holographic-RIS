%% 这个脚本准备尝试对多用户进行符号解调
% 这个脚本实现了对多用户的符号解调，效果不错

clc; clear; close all;

%% ================= 多用户发送端参数设置 =================

U = 3;                       % 用户数
N_sym = 50;                 % 每个用户的符号数
M_mod = 4;                  % QPSK
bits_per_sym = log2(M_mod);

Rs = 1e6;                    % 符号率
sps = 8;                     % 每符号基带采样点
rolloff = 0.25;              % RRC滚降系数
span = 6;                    % RRC滤波器跨度，单位：符号

interp = 20;                 % 额外插值倍数
fc = 20e6;                   % 仿真载波频率

fs_bb = Rs * sps;            % 基带采样率
fs = fs_bb * interp;         % 插值后的采样率

rng(64);

%% ================= 多用户发送端主流程 =================

% 1. 生成多用户随机比特
bits_users = generate_multiuser_bits(U, N_sym, bits_per_sym);

% 2. 每个用户做 QPSK 映射
symbols_users = qpsk_modulate_multiuser_bits(bits_users);

% 3. 生成 RRC 成形滤波器
rrc = generate_rrc_filter(rolloff, span, sps);

% 4. 每个用户做 RRC 成形，得到复基带波形
tx_bb_users = pulse_shape_multiuser_symbols(symbols_users, rrc, sps);

% 5. 每个用户插值到高采样率
tx_bb_hi_users = interpolate_multiuser_baseband(tx_bb_users, interp);

snr = 10;
tx_bb_hi_users = awgn(tx_bb_hi_users,snr);

% 6. 构造时间轴
t = generate_time_axis(size(tx_bb_hi_users, 1), fs);

% 7. 每个用户上变频到实射频信号
tx_rf_users = upconvert_multiuser_to_rf(tx_bb_hi_users, fc, fs);

%% ================= 输出基本信息 =================

fprintf('Number of users          = %d\n', U);
fprintf('Symbols per user         = %d\n', N_sym);
fprintf('Bits per user            = %d\n', N_sym * bits_per_sym);
fprintf('Symbol rate              = %.2f MHz\n', Rs/1e6);
fprintf('Baseband sample rate     = %.2f MHz\n', fs_bb/1e6);
fprintf('RF simulation rate       = %.2f MHz\n', fs/1e6);
fprintf('Carrier frequency        = %.2f MHz\n', fc/1e6);
fprintf('Length of tx_bb_hi_users = %d samples\n', size(tx_bb_hi_users, 1));

%% ================= 多用户发送端可视化 =================

plot_multiuser_transmitter_waveforms( ...
    t, tx_bb_hi_users, tx_rf_users, symbols_users, fc);


%% ================= 多用户接收端参数设置 =================

RIS_row = 8;
RIS_col = 8;
M_ris = RIS_row * RIS_col;

lambda = 1;
d = 0.5 * lambda;

% 三个用户的真实到达角
% 每一行对应一个用户：[theta, phi]
user_angles_deg = [
     20,  10;
    -15,   5;
     35, -12
];

theta_users = user_angles_deg(:,1) * pi/180;
phi_users   = user_angles_deg(:,2) * pi/180;

ref_amp = 1.5;

lpf_cutoff = 5e6;
fir_order = 800;

constellation = qammod((0:M_mod-1).', M_mod, 'gray', ...
    'UnitAveragePower', true);

%% ================= 多用户接收端主流程 =================

% 1. 生成多用户到 RIS 阵列的阵列流形矩阵
% H_ris 的尺寸为 M_ris x U
H_ris = generate_multiuser_ris_steering_matrix( ...
    RIS_row, RIS_col, theta_users, phi_users, d, lambda);

% 2. 生成每个 RIS 单元的参考符号
% 参考信号是接收端本地已知的，不区分用户
ref_symbols_ris = generate_reference_symbols_ris( ...
    N_sym, M_ris, ref_amp);

% 3. 计算符号能量积分窗口
[sym_center_idx, win_len] = get_symbol_windows( ...
    N_sym, sps, interp, span);

% 4. 仿真多用户 SIMO 符号级能量接收机
E_simo = simulate_multiuser_simo_symbol_energy_receiver( ...
    tx_bb_hi_users, H_ris, ref_symbols_ris, rrc, sps, interp, ...
    fc, fs, sym_center_idx, win_len, lpf_cutoff, fir_order);

% 5. 估计符号窗口内的脉冲能量系数
Cp = estimate_pulse_energy_coefficient( ...
    N_sym, rrc, sps, interp, sym_center_idx, win_len);

%% ================= 多用户检测算法对比：ML / GS / GN =================

% ---------- 1. 多用户联合 ML 检测 ----------
symbols_users_est_ml = multiuser_symbol_energy_ml_detector( ...
    E_simo, H_ris, ref_symbols_ris, constellation, Cp);

% ---------- 2. 多用户逐符号 GS 检测 ----------
t0_gs = 1000;

symbols_users_est_gs = multiuser_symbol_energy_gs_detector( ...
    E_simo, H_ris, ref_symbols_ris, constellation, Cp, t0_gs);

% ---------- 3. 多用户逐符号 GN 检测 ----------
maxIter_gn = 200;

[symbols_users_est_gn, symbols_users_cont_gn] = multiuser_symbol_energy_gn_detector( ...
    E_simo, H_ris, ref_symbols_ris, constellation, Cp, maxIter_gn, N_sym, U);

% ---------- 4. 计算每个用户的 SER ----------
ser_users_ml = mean(symbols_users_est_ml ~= symbols_users, 1);
ser_users_gs = mean(symbols_users_est_gs ~= symbols_users, 1);
ser_users_gn = mean(symbols_users_est_gn ~= symbols_users, 1);

fprintf('\n========== Multiuser Detection Comparison ==========\n');

for u = 1:U
    fprintf('User %d ML SER = %.4f, GS SER = %.4f, GN SER = %.4f\n', ...
        u, ser_users_ml(u), ser_users_gs(u), ser_users_gn(u));
end

fprintf('Average ML SER = %.4f\n', mean(ser_users_ml));
fprintf('Average GS SER = %.4f\n', mean(ser_users_gs));
fprintf('Average GN SER = %.4f\n', mean(ser_users_gn));

%% ================= 结果可视化 =================

figure;
for u = 1:U
    subplot(1,U,u);
    plot(real(symbols_users(:,u)), imag(symbols_users(:,u)), ...
        'ko', 'LineWidth', 1.5); hold on;
    plot(real(symbols_users_cont_gn(:,u)), imag(symbols_users_cont_gn(:,u)), ...
        'rx', 'LineWidth', 1.2);
    grid on; axis equal;
    xlabel('In-phase');
    ylabel('Quadrature');
    title(sprintf('User %d GN Continuous Estimates', u));
    legend('True symbols', 'GN continuous');
end

%% ================= 检测结果可视化 =================

figure;

for u = 1:U
    subplot(U, 1, u);

    plot(real(symbols_users(:,u)), imag(symbols_users(:,u)), ...
        'ko', 'LineWidth', 1.5); hold on;

    plot(real(symbols_users_est_ml(:,u)), imag(symbols_users_est_ml(:,u)), ...
        'rx', 'LineWidth', 1.2);

    plot(real(symbols_users_est_gs(:,u)), imag(symbols_users_est_gs(:,u)), ...
        'b+', 'LineWidth', 1.2);

    plot(real(symbols_users_est_gn(:,u)), imag(symbols_users_est_gn(:,u)), ...
        'ms', 'LineWidth', 1.2);

    grid on; axis equal;
    xlabel('In-phase');
    ylabel('Quadrature');
    title(sprintf('User %d Detection Comparison', u));
    legend('True', 'ML', 'GS', 'GN');
end

figure;
imagesc(E_simo);
colorbar;
xlabel('Symbol index');
ylabel('RIS element index');
title('Multiuser SIMO Symbol Energy Matrix E_{m,k}');

%% ========================================================================
%                         多用户发送端函数区
% ========================================================================

function bits_users = generate_multiuser_bits(U, N_sym, bits_per_sym)
%GENERATE_MULTIUSER_BITS 生成多用户随机比特
%
% 输入：
%   U            : 用户数
%   N_sym        : 每个用户符号数
%   bits_per_sym : 每个符号对应的比特数
%
% 输出：
%   bits_users   : N_bits x U，比特矩阵，每一列对应一个用户

    N_bits = N_sym * bits_per_sym;

    bits_users = randi([0 1], N_bits, U);
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
%   span    : 滤波器跨度，单位为符号
%   sps     : 每符号采样点数
%
% 输出：
%   rrc     : RRC滤波器系数

    rrc = rcosdesign(rolloff, span, sps, 'sqrt');

    rrc = rrc(:);
end


function tx_bb_users = pulse_shape_multiuser_symbols(symbols_users, rrc, sps)
%PULSE_SHAPE_MULTIUSER_SYMBOLS 对每个用户的符号做 RRC 成形
%
% 输入：
%   symbols_users : N_sym x U，多用户符号矩阵
%   rrc           : RRC 滤波器
%   sps           : 每符号采样点数
%
% 输出：
%   tx_bb_users   : N_sample x U，每一列为一个用户的复基带波形

    [~, U] = size(symbols_users);

    % 先处理第一个用户，确定输出长度
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
%INTERPOLATE_MULTIUSER_BASEBAND 对多用户基带波形插值
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
%   t        : 时间列向量

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


function plot_multiuser_transmitter_waveforms( ...
    t, tx_bb_hi_users, tx_rf_users, symbols_users, fc)
%PLOT_MULTIUSER_TRANSMITTER_WAVEFORMS 绘制多用户发送端关键波形

    [~, U] = size(tx_bb_hi_users);

    figure('Position', [100, 100, 1200, 850]);

    subplot(4,1,1);
    hold on;
    for u = 1:U
        plot(real(symbols_users(:,u)), imag(symbols_users(:,u)), ...
            'o', 'LineWidth', 1.2);
    end
    grid on; axis equal;
    xlabel('In-phase');
    ylabel('Quadrature');
    title('Multiuser QPSK Constellation Symbols');
    legend(compose('User %d', 1:U));

    subplot(4,1,2);
    hold on;
    for u = 1:U
        plot(t*1e6, real(tx_bb_hi_users(:,u)), 'LineWidth', 1.0);
    end
    grid on;
    xlabel('Time / \mus');
    ylabel('Amplitude');
    title('Real Part of Pulse-shaped Baseband Signals');
    legend(compose('User %d', 1:U));

    subplot(4,1,3);
    hold on;
    for u = 1:U
        plot(t*1e6, abs(tx_bb_hi_users(:,u)), 'LineWidth', 1.0);
    end
    grid on;
    xlabel('Time / \mus');
    ylabel('Envelope');
    title('Baseband Envelopes |x_{BB,u}(t)|');
    legend(compose('User %d', 1:U));

    subplot(4,1,4);

    if length(t) > 1
        T_show = max(2e-6, t(end));
        N_show = find(t <= T_show, 1, 'last');
        if isempty(N_show)
%             N_show = min(length(t), 2000);
            N_show = length(t);
        end
    else
        N_show = 1;
    end

    hold on;
    for u = 1:U
        plot(t(1:N_show)*1e6, tx_rf_users(1:N_show,u), ...
            'LineWidth', 1.0);
    end
    grid on;
    xlabel('Time / \mus');
    ylabel('Amplitude');
    title('Upconverted Real RF Signals');
    legend(compose('User %d', 1:U));

    sgtitle('Multiuser Transmitter Signal Generation');

    %% ================= 多用户发送波形叠加 =================
    
    tx_bb_hi_sum = sum(tx_bb_hi_users, 2);   % 多用户复基带叠加波形
    tx_rf_sum    = sum(tx_rf_users, 2);      % 多用户实射频叠加波形

    figure;

    subplot(3,1,1);
    plot(t*1e6, real(tx_bb_hi_sum), 'LineWidth', 1.2); hold on;
    plot(t*1e6, imag(tx_bb_hi_sum), 'LineWidth', 1.2);
    grid on;
    xlabel('Time / \mus');
    ylabel('Amplitude');
    legend('I', 'Q');
    title('Superposed Multiuser Complex Baseband Signal');
    
    subplot(3,1,2);
    plot(t*1e6, abs(tx_bb_hi_sum), 'LineWidth', 1.2);
    grid on;
    xlabel('Time / \mus');
    ylabel('Envelope');
    title('Envelope of Superposed Multiuser Baseband Signal');
    
    subplot(3,1,3);
    N_show = max(length(t), 2000);
    plot(t(1:N_show)*1e6, tx_rf_sum(1:N_show), 'LineWidth', 1.0);
    grid on;
    xlabel('Time / \mus');
    ylabel('Amplitude');
    title('Superposed Multiuser Real RF Signal');

end


%% ========================================================================
%                         多用户接收端函数区
% ========================================================================

function H_ris = generate_multiuser_ris_steering_matrix( ...
    RIS_row, RIS_col, theta_users, phi_users, d, lambda)
%GENERATE_MULTIUSER_RIS_STEERING_MATRIX
% 生成多用户到 RIS 阵列的阵列流形矩阵。
%
% 输入：
%   RIS_row     : RIS 行数
%   RIS_col     : RIS 列数
%   theta_users : U x 1，每个用户的 theta
%   phi_users   : U x 1，每个用户的 phi
%   d           : 阵元间距
%   lambda      : 波长
%
% 输出：
%   H_ris       : M_ris x U，H_ris(m,u) 为第 u 个用户到第 m 个 RIS 单元的阵列系数

    theta_users = theta_users(:);
    phi_users = phi_users(:);

    U = length(theta_users);

    if length(phi_users) ~= U
        error('theta_users 和 phi_users 的长度必须一致。');
    end

    M_ris = RIS_row * RIS_col;

    H_ris = zeros(M_ris, U);

    for u = 1:U
        H_ris(:,u) = generate_ris_steering_vector( ...
            RIS_row, RIS_col, theta_users(u), phi_users(u), d, lambda);
    end
end


function h_ris = generate_ris_steering_vector(RIS_row, RIS_col, theta, phi, d, lambda)
%GENERATE_RIS_STEERING_VECTOR
% 生成二维 RIS 平面阵列的阵列流形。
%
% 假设 RIS 位于 y-z 平面，x 轴为法向。
% row 对应 z 方向，col 对应 y 方向。
%
% 输入：
%   RIS_row : RIS 行数
%   RIS_col : RIS 列数
%   theta   : 入射角 theta，单位 rad
%   phi     : 入射角 phi，单位 rad
%   d       : 阵元间距
%   lambda  : 波长
%
% 输出：
%   h_ris   : M_ris x 1 的复阵列流形向量

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
%GENERATE_REFERENCE_SYMBOLS_RIS
% 为每个 RIS 单元生成已知参考符号序列。
%
% 这里采用：
%   公共确定性相位序列 + 每个 RIS 单元独立相位偏置。
%
% 输入：
%   N_sym   : 符号数
%   M_ris   : RIS 单元数
%   ref_amp : 参考信号幅度
%
% 输出：
%   ref_symbols_ris : N_sym x M_ris 的参考符号矩阵

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
%GET_SYMBOL_WINDOWS
% 计算每个符号的中心采样点和积分窗口长度。
%
% 输入：
%   N_sym  : 符号数
%   sps    : 每符号基带采样点数
%   interp : 插值倍数
%   span   : RRC 滤波器跨度
%
% 输出：
%   sym_center_idx : N_sym x 1，每个符号中心索引
%   win_len        : 每个符号能量积分窗口长度

    sym_samp_hi = sps * interp;

    rrc_delay_bb = span * sps / 2;
    rrc_delay_hi = rrc_delay_bb * interp;

    sym_center_idx = round(rrc_delay_hi + 1 + (0:N_sym-1).' * sym_samp_hi);

    win_len = sym_samp_hi;
end


function E_simo = simulate_multiuser_simo_symbol_energy_receiver( ...
    tx_bb_hi_users, H_ris, ref_symbols_ris, rrc, sps, interp, ...
    fc, fs, sym_center_idx, win_len, lpf_cutoff, fir_order)
%SIMULATE_MULTIUSER_SIMO_SYMBOL_ENERGY_RECEIVER
% 仿真多用户 SIMO-RIS 符号级能量接收机。
%
% 对每个 RIS 单元 m：
%   1. 多个用户经阵列流形叠加：
%        x_m(t) = sum_u h_{m,u} x_u(t)
%   2. 与第 m 个单元的本地参考信号叠加
%   3. 上变频为实射频信号
%   4. 平方律检波 + 低通
%   5. 每个符号周期提取平均能量
%
% 输入：
%   tx_bb_hi_users : N_sample x U，多用户复基带波形
%   H_ris          : M_ris x U，多用户阵列流形矩阵
%   ref_symbols_ris: N_sym x M_ris，RIS 参考符号
%   rrc, sps, interp : 参考信号成形参数
%   fc, fs         : 载波频率和采样率
%   sym_center_idx : N_sym x 1，符号中心索引
%   win_len        : 符号积分窗口长度
%   lpf_cutoff     : 平方律检波低通截止频率
%   fir_order      : 低通滤波器阶数
%
% 输出：
%   E_simo         : M_ris x N_sym，符号级能量矩阵

    [N_sample, U] = size(tx_bb_hi_users);
    [M_ris, U_h] = size(H_ris);
    [N_sym, M_ref] = size(ref_symbols_ris);

    if U_h ~= U
        error('H_ris 的列数必须等于用户数 U。');
    end

    if M_ref ~= M_ris
        error('ref_symbols_ris 的列数必须等于 RIS 单元数。');
    end

    E_simo = zeros(M_ris, N_sym);

    t = (0:N_sample-1).' / fs;
    carrier = exp(1j * 2*pi*fc*t);

    for m = 1:M_ris

        % 多用户在第 m 个 RIS 单元处的未知复基带叠加信号
        rx_unknown_bb_m = zeros(N_sample, 1);

        for u = 1:U
            rx_unknown_bb_m = rx_unknown_bb_m ...
                + H_ris(m,u) * tx_bb_hi_users(:,u);
        end

        % 第 m 个 RIS 单元的参考复基带信号
        ref_bb_m = pulse_shape_symbols(ref_symbols_ris(:,m), rrc, sps);
        ref_bb_hi_m = interpolate_baseband(ref_bb_m, interp);

        L = min([length(rx_unknown_bb_m), length(ref_bb_hi_m), length(carrier)]);

        rx_unknown_bb_m = rx_unknown_bb_m(1:L);
        ref_bb_hi_m = ref_bb_hi_m(1:L);
        carrier_m = carrier(1:L);

        % 分别上变频为实射频信号
        rx_unknown_rf_m = real(rx_unknown_bb_m .* carrier_m);
        ref_rf_m = real(ref_bb_hi_m .* carrier_m);

        % 射频叠加
        mixed_rf_m = rx_unknown_rf_m + ref_rf_m;

        % 平方律检波 + 低通
        energy_waveform_m = square_law_detector( ...
            mixed_rf_m, fs, lpf_cutoff, fir_order);

        % 每个符号周期提取平均能量
        E_simo(m,:) = extract_symbol_energy_from_waveform( ...
            energy_waveform_m, sym_center_idx, win_len);
    end
end


function symbols_users_est = multiuser_symbol_energy_ml_detector( ...
    E_simo, H_ris, ref_symbols_ris, constellation, Cp)
%MULTIUSER_SYMBOL_ENERGY_ML_DETECTOR
% 多用户符号级能量联合 ML 检测器。
%
% 简化符号级模型：
%   E_{m,k} ≈ Cp * | sum_u h_{m,u} c_u + b_{m,k} |^2
%
% 对每个符号 k，遍历所有多用户星座组合：
%   [c_1, ..., c_U] ∈ C^U
%
% 选择预测能量与观测能量最接近的组合。
%
% 输入：
%   E_simo          : M_ris x N_sym，符号能量观测
%   H_ris           : M_ris x U，多用户阵列流形矩阵
%   ref_symbols_ris : N_sym x M_ris，参考符号矩阵
%   constellation   : 星座点列向量
%   Cp              : 脉冲能量系数
%
% 输出：
%   symbols_users_est : N_sym x U，估计的多用户符号

    constellation = constellation(:);

    [M_ris, N_sym] = size(E_simo);
    [M_h, U] = size(H_ris);

    if M_h ~= M_ris
        error('H_ris 的行数必须等于 E_simo 的行数。');
    end

    if size(ref_symbols_ris,1) ~= N_sym || size(ref_symbols_ris,2) ~= M_ris
        error('ref_symbols_ris 的尺寸必须为 N_sym x M_ris。');
    end

    % 生成所有多用户星座组合
    combo_symbols = generate_constellation_combinations(constellation, U);
    N_combo = size(combo_symbols, 1);

    symbols_users_est = zeros(N_sym, U);

    for k = 1:N_sym

        metric = zeros(N_combo, 1);

        b_k = ref_symbols_ris(k,:).';   % M_ris x 1

        for cc = 1:N_combo

            c_vec = combo_symbols(cc,:).';   % U x 1

            % 多用户在每个 RIS 单元上的叠加符号
            y_pred = H_ris * c_vec + b_k;    % M_ris x 1

            E_pred = Cp * abs(y_pred).^2;

            metric(cc) = sum((E_simo(:,k) - E_pred).^2);
        end

        [~, id_min] = min(metric);

        symbols_users_est(k,:) = combo_symbols(id_min,:);
    end
end


function combo_symbols = generate_constellation_combinations(constellation, U)
%GENERATE_CONSTELLATION_COMBINATIONS
% 生成 U 个用户的所有星座组合。
%
% 输入：
%   constellation : 星座点列向量，长度为 M_mod
%   U             : 用户数
%
% 输出：
%   combo_symbols : M_mod^U x U，每一行是一个多用户符号组合

    constellation = constellation(:);
    M_mod = length(constellation);

    grid_cell = cell(1, U);

    [grid_cell{:}] = ndgrid(1:M_mod);

    N_combo = M_mod^U;

    combo_symbols = zeros(N_combo, U);

    for u = 1:U
        idx_u = grid_cell{u}(:);
        combo_symbols(:,u) = constellation(idx_u);
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


function energy_est = square_law_detector(rx_rf, fs, lpf_cutoff, fir_order)
%SQUARE_LAW_DETECTOR
% 平方律检波器。
%
% 对实射频信号：
%   rx_rf(t) = Re{u(t) exp(j2πf_ct)}
%
% 平方后低通：
%   LPF{rx_rf^2} ≈ 1/2 |u(t)|^2
%
% 因此最后乘以 2 得到基带能量估计。

    rx_rf = rx_rf(:);

    rx_square = rx_rf.^2;

    b_lpf = fir1(fir_order, lpf_cutoff/(fs/2), 'low');

    rx_energy_lpf = filtfilt(b_lpf, 1, rx_square);

    energy_est = 2 * rx_energy_lpf;

    energy_est = max(real(energy_est), 0);

    energy_est = energy_est(:);
end


function E_symbol = extract_symbol_energy_from_waveform(energy_waveform, sym_center_idx, win_len)
%EXTRACT_SYMBOL_ENERGY_FROM_WAVEFORM
% 从能量波形中提取每个符号平均能量。

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
%PULSE_SHAPE_SYMBOLS
% 对符号序列进行上采样和 RRC 成形滤波。
%
% 输入：
%   symbols : 复符号列向量
%   rrc     : RRC 滤波器
%   sps     : 每符号采样点数
%
% 输出：
%   tx_bb   : 成形后的复基带波形

    symbols = symbols(:);

    tx_bb = upfirdn(symbols, rrc, sps, 1);

    tx_bb = tx_bb(:);
end


function tx_bb_hi = interpolate_baseband(tx_bb, interp)
%INTERPOLATE_BASEBAND
% 对基带信号进行插值。

    tx_bb = tx_bb(:);

    if interp == 1
        tx_bb_hi = tx_bb;
    else
        tx_bb_hi = resample(tx_bb, interp, 1);
    end

    tx_bb_hi = tx_bb_hi(:);
end

function symbols_users_est = multiuser_symbol_energy_gs_detector( ...
    E_simo, H_ris, ref_symbols_ris, constellation, Cp, t0)
%MULTIUSER_SYMBOL_ENERGY_GS_DETECTOR
% 使用旧 biased_gs_algorithm 对多用户符号级能量模型做逐符号恢复。
%
% 符号级能量模型：
%   E_{m,k} ≈ Cp * | sum_u h_{m,u} s_{u,k} + b_{m,k} |^2
%
% 转换为幅值模型：
%   sqrt(E_{m,k}) ≈ | sqrt(Cp) * H_ris * s_k
%                      + sqrt(Cp) * b_k |
%
% 如果 biased_gs_algorithm 的接口是：
%   z = abs(A' * s + b)
%
% 则应传入：
%   A_gs = (sqrt(Cp) * H_ris).'
%
% 输入：
%   E_simo          : M_ris x N_sym，符号级能量观测
%   H_ris           : M_ris x U，多用户阵列流形矩阵
%   ref_symbols_ris : N_sym x M_ris，参考符号矩阵
%   constellation   : 星座点列向量
%   Cp              : 符号窗口脉冲能量系数
%   t0              : GS 迭代次数
%
% 输出：
%   symbols_users_est : N_sym x U，估计的多用户符号

    constellation = constellation(:);

    [M_ris, N_sym] = size(E_simo);
    [M_h, U] = size(H_ris);

    if M_h ~= M_ris
        error('H_ris 的行数必须等于 E_simo 的行数。');
    end

    if size(ref_symbols_ris,1) ~= N_sym || size(ref_symbols_ris,2) ~= M_ris
        error('ref_symbols_ris 的尺寸必须为 N_sym x M_ris。');
    end

    if Cp <= 0
        error('Cp 必须为正数。');
    end

    symbols_users_est = zeros(N_sym, U);

    % 幅值模型中的多用户观测矩阵
    A_symbol = sqrt(Cp) * H_ris;    % M_ris x U

    % 旧 GS 函数使用 z = abs(A' * s + b)
    A_gs = A_symbol.';              % U x M_ris

    for k = 1:N_sym

        % 第 k 个符号时刻的幅值观测
        z_gs = sqrt(max(E_simo(:,k), 0));   % M_ris x 1

        % 第 k 个符号时刻的参考偏置
        b_gs = sqrt(Cp) * ref_symbols_ris(k,:).';   % M_ris x 1

        % GS 连续恢复 U 个用户的符号向量
        s_cont = biased_gs_algorithm(z_gs, A_gs, b_gs, t0);

        s_cont = s_cont(:);

        % 每个用户分别映射到最近星座点
        for u = 1:U
            [~, id_min] = min(abs(s_cont(u) - constellation));
            symbols_users_est(k,u) = constellation(id_min);
        end
    end
end

function [symbols_users_est, symbols_users_cont] = multiuser_symbol_energy_gn_detector( ...
    E_simo, H_ris, ref_symbols_ris, constellation, Cp, maxIter,N_sym,U)
%MULTIUSER_SYMBOL_ENERGY_GN_DETECTOR
% 使用 biased_gn_algorithm 对多用户符号级能量模型做逐符号恢复。
%
% 符号级能量模型：
%   E_{m,k} ≈ Cp * | sum_u h_{m,u} s_{u,k} + b_{m,k} |^2
%
% 转换为幅值模型：
%   sqrt(E_{m,k}) ≈ | sqrt(Cp) * H_ris * s_k
%                      + sqrt(Cp) * b_k |
%
% biased_gn_algorithm 的接口应为：
%   z = abs(A*s + b)
%
% 输入：
%   E_simo          : M_ris x N_sym，符号级能量观测
%   H_ris           : M_ris x U，多用户阵列流形矩阵
%   ref_symbols_ris : N_sym x M_ris，参考符号矩阵
%   constellation   : 星座点列向量
%   Cp              : 脉冲能量系数
%   maxIter         : GN 最大迭代次数
%
% 输出：
%   symbols_users_est : N_sym x U，估计的多用户符号

    symbols_users_est = zeros(N_sym, U);
    symbols_users_cont = zeros(N_sym, U);

    constellation = constellation(:);

    [M_ris, N_sym] = size(E_simo);
    [M_h, U] = size(H_ris);

    if M_h ~= M_ris
        error('H_ris 的行数必须等于 E_simo 的行数。');
    end

    if size(ref_symbols_ris,1) ~= N_sym || size(ref_symbols_ris,2) ~= M_ris
        error('ref_symbols_ris 的尺寸必须为 N_sym x M_ris。');
    end

    if Cp <= 0
        error('Cp 必须为正数。');
    end

    symbols_users_est = zeros(N_sym, U);

    % 幅值模型中的多用户观测矩阵
    A_gn = sqrt(Cp) * H_ris;     % M_ris x U

    for k = 1:N_sym

        % 第 k 个符号时刻的幅值观测
        z_gn = sqrt(max(E_simo(:,k), 0));   % M_ris x 1

        % 第 k 个符号时刻的参考偏置
        b_gn = sqrt(Cp) * ref_symbols_ris(k,:).';   % M_ris x 1

        % GN 连续恢复 U 个用户的符号向量
        s_cont = biased_gn_algorithm(z_gn, A_gn, b_gn, maxIter);

        s_cont = s_cont(:);

        % 每个用户分别映射到最近星座点
        for u = 1:U
            [~, id_min] = min(abs(s_cont(u) - constellation));
            symbols_users_est(k,u) = constellation(id_min);
        end

        symbols_users_cont(k,:) = s_cont(:).';
    end
end