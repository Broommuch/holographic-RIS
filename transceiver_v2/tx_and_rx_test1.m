clc; clear; close all;

%% ================= 发送端参数设置 =================

N_sym = 50;              % 发送符号数
M_mod = 4;                % QPSK
bits_per_sym = log2(M_mod);

Rs = 1e6;                 % 符号率
sps = 8;                  % 每符号基带采样点
rolloff = 0.25;           % RRC滚降系数
span = 6;                 % RRC滤波器跨度，单位：符号

interp = 20;              % 额外插值倍数，用于生成更平滑的射频波形
fc = 20e6;                % 仿真载波频率

fs_bb = Rs * sps;         % 基带采样率
fs = fs_bb * interp;      % 插值后的采样率

rng(64);

%% ================= 发送端主流程 =================

% 1. 生成随机比特
bits_tx = generate_random_bits(N_sym, bits_per_sym);

% 2. QPSK星座映射
symbols_tx = qpsk_modulate_bits(bits_tx);

% 3. 生成RRC滤波器
rrc = generate_rrc_filter(rolloff, span, sps);

% 4. RRC成形滤波，得到基带波形
tx_bb = pulse_shape_symbols(symbols_tx, rrc, sps);

% 5. 插值，提高射频仿真采样率
tx_bb_hi = interpolate_baseband(tx_bb, interp);
snr = 10;
tx_bb_hi = awgn(tx_bb_hi,snr);

% 6. 构造时间轴
t = generate_time_axis(length(tx_bb_hi), fs);

% 7. 上变频到射频
tx_rf = upconvert_to_rf(tx_bb_hi, fc, fs);

%% ================= 发送端信息输出 =================

fprintf('Number of symbols       = %d\n', N_sym);
fprintf('Number of bits          = %d\n', length(bits_tx));
fprintf('Symbol rate             = %.2f MHz\n', Rs/1e6);
fprintf('Baseband sample rate    = %.2f MHz\n', fs_bb/1e6);
fprintf('RF simulation rate      = %.2f MHz\n', fs/1e6);
fprintf('Carrier frequency       = %.2f MHz\n', fc/1e6);
fprintf('Length of tx_bb         = %d samples\n', length(tx_bb));
fprintf('Length of tx_bb_hi      = %d samples\n', length(tx_bb_hi));
fprintf('Length of tx_rf         = %d samples\n', length(tx_rf));

%% ================= 发送端可视化 =================

plot_transmitter_waveforms(t, tx_bb_hi, tx_rf, symbols_tx, fc);


%% ================= 接收端参数设置 =================

RIS_row = 8;
RIS_col = 8;
M_ris = RIS_row * RIS_col;

lambda = 1;
d = 0.5 * lambda;

theta_u = 20 * pi/180;      % 用户到达角，按你的坐标系修改
phi_u   = 10 * pi/180;

ref_amp = 1.5;              % 参考信号幅度

lpf_cutoff = 5e6;           % 平方律检波后的低通截止频率
fir_order = 800;            % 低通滤波器阶数

M_mod = 4;
constellation = qammod((0:M_mod-1).', M_mod, 'gray', ...
    'UnitAveragePower', true);

%% ================= 接收端主流程 =================

% 1. 生成 RIS 阵列流形
h_ris = generate_ris_steering_vector( ...
    RIS_row, RIS_col, theta_u, phi_u, d, lambda);

% 2. 为每个 RIS 单元生成符号级参考序列
ref_symbols_ris = generate_reference_symbols_ris( ...
    N_sym, M_ris, ref_amp);

% 3. 估计符号窗口中心和窗口长度
[sym_center_idx, win_len] = get_symbol_windows( ...
    N_sym, sps, interp, span);

% 4. 仿真每个 RIS 单元的接收、参考叠加、平方律检波、符号能量提取
E_simo = simulate_simo_symbol_energy_receiver( ...
    tx_bb_hi, h_ris, ref_symbols_ris, rrc, sps, interp, ...
    fc, fs, sym_center_idx, win_len, lpf_cutoff, fir_order);

% 5. 估计成形脉冲在符号窗口内的能量系数
Cp = estimate_pulse_energy_coefficient( ...
    N_sym, rrc, sps, interp, sym_center_idx, win_len);

% 6. 符号级 SIMO 能量 ML 检测
symbols_est = symbol_energy_ml_detector( ...
    E_simo, h_ris, ref_symbols_ris, constellation, Cp);

% 7. 计算符号错误率
ser = mean(symbols_est(:) ~= symbols_tx(:));

fprintf('Symbol-level SIMO energy SER = %.4f\n', ser);

%% ================= 接收端结果可视化 =================

figure;
plot(real(symbols_tx), imag(symbols_tx), 'o', 'LineWidth', 1.5); hold on;
plot(real(symbols_est), imag(symbols_est), 'x', 'LineWidth', 1.5);
grid on; axis equal;
xlabel('In-phase');
ylabel('Quadrature');
legend('True QPSK symbols', 'Estimated QPSK symbols');
title('Symbol-level SIMO Energy Detection Result');

figure;
imagesc(E_simo);
colorbar;
xlabel('Symbol index');
ylabel('RIS element index');
title('Extracted Symbol Energy Matrix E_{m,k}');


%% ========================================================================
%                              函数区
% ========================================================================

function bits = generate_random_bits(N_sym, bits_per_sym)
%GENERATE_RANDOM_BITS 生成随机比特序列
%
% 输入：
%   N_sym        : 符号数
%   bits_per_sym : 每个符号对应的比特数
%
% 输出：
%   bits         : 随机比特列向量

    N_bits = N_sym * bits_per_sym;
    bits = randi([0 1], N_bits, 1);
end


function symbols = qpsk_modulate_bits(bits)
%QPSK_MODULATE_BITS 将比特映射为QPSK星座符号
%
% 使用Gray映射，单位平均功率归一化。
%
% 输入：
%   bits    : 比特列向量，长度应为偶数
%
% 输出：
%   symbols : QPSK复符号列向量

    bits = bits(:);

    if mod(length(bits), 2) ~= 0
        error('QPSK调制要求比特数为2的整数倍。');
    end

    bits_reshape = reshape(bits, 2, []).';

    symbol_idx = bi2de(bits_reshape, 'left-msb');

    symbols = qammod(symbol_idx, 4, 'gray', ...
        'UnitAveragePower', true);

    symbols = symbols(:);
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


function tx_bb = pulse_shape_symbols(symbols, rrc, sps)
%PULSE_SHAPE_SYMBOLS 对符号序列进行上采样和RRC成形滤波
%
% 输入：
%   symbols : 复符号列向量
%   rrc     : RRC滤波器
%   sps     : 每符号采样点数
%
% 输出：
%   tx_bb   : 成形后的复基带波形

    symbols = symbols(:);

    tx_bb = upfirdn(symbols, rrc, sps, 1);

    tx_bb = tx_bb(:);
end


function tx_bb_hi = interpolate_baseband(tx_bb, interp)
%INTERPOLATE_BASEBAND 对基带信号进行插值
%
% 输入：
%   tx_bb    : 原始成形基带波形
%   interp   : 插值倍数
%
% 输出：
%   tx_bb_hi : 插值后的高采样率基带波形

    tx_bb = tx_bb(:);

    if interp == 1
        tx_bb_hi = tx_bb;
    else
        tx_bb_hi = resample(tx_bb, interp, 1);
    end

    tx_bb_hi = tx_bb_hi(:);
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


function tx_rf = upconvert_to_rf(tx_bb, fc, fs)
%UPCONVERT_TO_RF 将复基带信号上变频为实射频信号
%
% 复基带模型：
%   x_RF(t) = Re{ x_BB(t) exp(j2πf_ct) }
%
% 输入：
%   tx_bb : 复基带波形
%   fc    : 载波频率
%   fs    : 采样率
%
% 输出：
%   tx_rf : 实射频波形

    tx_bb = tx_bb(:);

    t = (0:length(tx_bb)-1).' / fs;

    carrier = exp(1j * 2*pi*fc*t);

    tx_rf = real(tx_bb .* carrier);

    tx_rf = tx_rf(:);
end


function plot_transmitter_waveforms(t, tx_bb_hi, tx_rf, symbols_tx, fc)
%PLOT_TRANSMITTER_WAVEFORMS 绘制发送端关键波形

    figure('Position', [100, 100, 1200, 850]);

    subplot(4,1,1);
    plot(real(symbols_tx), imag(symbols_tx), 'o', 'LineWidth', 1.5);
    grid on; axis equal;
    xlabel('In-phase');
    ylabel('Quadrature');
    title('QPSK Constellation Symbols');

    subplot(4,1,2);
    plot(t*1e6, real(tx_bb_hi), 'b'); hold on;
    plot(t*1e6, imag(tx_bb_hi), 'r');
    grid on;
    xlabel('Time / \mus');
    ylabel('Amplitude');
    legend('I', 'Q');
    title('Pulse-shaped Complex Baseband Signal');

    subplot(4,1,3);
    plot(t*1e6, abs(tx_bb_hi), 'k', 'LineWidth', 1.2);
    grid on;
    xlabel('Time / \mus');
    ylabel('Envelope');
    title('Baseband Envelope |x_{BB}(t)|');

    subplot(4,1,4);
    N_show = min(length(t), round(5 / fc * length(t) / (t(end)-t(1))));
    if N_show < 100
        N_show = min(length(t), 2000);
    end

    plot(t(1:N_show)*1e6, tx_rf(1:N_show), 'LineWidth', 1.0);
    grid on;
    xlabel('Time / \mus');
    ylabel('Amplitude');
    title('Upconverted Real RF Signal');

    sgtitle('Transmitter Signal Generation');
end

%% 接收端函数
%% ========================================================================
%                         接收端函数区
% ========================================================================

function h_ris = generate_ris_steering_vector(RIS_row, RIS_col, theta, phi, d, lambda)
%GENERATE_RIS_STEERING_VECTOR 生成二维 RIS 平面阵列的阵列流形
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
%GENERATE_REFERENCE_SYMBOLS_RIS 为每个 RIS 单元生成已知参考符号序列
%
% 这里采用：
%   公共确定性相位序列 + 每个 RIS 单元独立相位偏置
%
% 输入：
%   N_sym   : 符号数
%   M_ris   : RIS 单元数
%   ref_amp : 参考信号幅度
%
% 输出：
%   ref_symbols_ris : N_sym x M_ris 的参考符号矩阵

    k = (0:N_sym-1).';

    % 公共参考相位序列，避免所有符号参考相位完全相同
    base_phase = mod(2*pi*0.137*k.^2 + pi/7*k, 2*pi);

    % 每个 RIS 单元一个独立参考相位偏置
    rng(2028);
    ris_phase = 2*pi*rand(1, M_ris);

    ref_symbols_ris = zeros(N_sym, M_ris);

    for m = 1:M_ris
        ref_symbols_ris(:,m) = ref_amp * exp(1j * (base_phase + ris_phase(m)));
    end
end


function [sym_center_idx, win_len] = get_symbol_windows(N_sym, sps, interp, span)
%GET_SYMBOL_WINDOWS 计算每个符号的中心采样点和积分窗口长度
%
% 输入：
%   N_sym  : 符号数
%   sps    : 每符号基带采样点数
%   interp : 插值倍数
%   span   : RRC 滤波器跨度
%
% 输出：
%   sym_center_idx : N_sym x 1，每个符号中心在高采样率波形中的索引
%   win_len        : 每个符号积分窗口长度

    sym_samp_hi = sps * interp;

    rrc_delay_bb = span * sps / 2;
    rrc_delay_hi = rrc_delay_bb * interp;

    sym_center_idx = round(rrc_delay_hi + 1 + (0:N_sym-1).' * sym_samp_hi);

    win_len = sym_samp_hi;
end


function E_simo = simulate_simo_symbol_energy_receiver( ...
    tx_bb_hi, h_ris, ref_symbols_ris, rrc, sps, interp, ...
    fc, fs, sym_center_idx, win_len, lpf_cutoff, fir_order)
%SIMULATE_SIMO_SYMBOL_ENERGY_RECEIVER
% 仿真 SIMO-RIS 符号级能量接收端。
%
% 对每个 RIS 单元：
%   1. 未知信号乘阵列流形 h_m
%   2. 生成该单元参考信号
%   3. 未知信号和参考信号在射频叠加
%   4. 平方律检波 + 低通滤波
%   5. 每个符号周期提取平均能量
%
% 输入：
%   tx_bb_hi        : 高采样率未知信号复基带波形
%   h_ris           : M_ris x 1 阵列流形
%   ref_symbols_ris : N_sym x M_ris 参考符号矩阵
%   rrc             : RRC 滤波器
%   sps             : 每符号基带采样点数
%   interp          : 插值倍数
%   fc              : 载波频率
%   fs              : 高采样率
%   sym_center_idx  : 符号中心索引
%   win_len         : 符号能量窗口长度
%   lpf_cutoff      : 平方律检波低通截止频率
%   fir_order       : 低通 FIR 阶数
%
% 输出：
%   E_simo          : M_ris x N_sym，每个 RIS 单元每个符号的平均能量

    tx_bb_hi = tx_bb_hi(:);
    h_ris = h_ris(:);

    [N_sym, M_ris] = size(ref_symbols_ris);

    E_simo = zeros(M_ris, N_sym);

    % 时间轴和载波
    t = (0:length(tx_bb_hi)-1).' / fs;
    carrier = exp(1j * 2*pi*fc*t);

    for m = 1:M_ris

        % 第 m 个 RIS 单元接收到的未知复基带信号
        rx_unknown_bb_m = h_ris(m) * tx_bb_hi;

        % 第 m 个 RIS 单元的参考复基带信号
        ref_bb_m = pulse_shape_symbols(ref_symbols_ris(:,m), rrc, sps);
        ref_bb_hi_m = interpolate_baseband(ref_bb_m, interp);

        % 对齐长度
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


function energy_est = square_law_detector(rx_rf, fs, lpf_cutoff, fir_order)
%SQUARE_LAW_DETECTOR 平方律检波器
%
% 对实射频信号：
%   rx_rf(t) = Re{u(t) exp(j2πf_ct)}
%
% 平方后低通：
%   LPF{rx_rf^2} ≈ 1/2 |u(t)|^2
%
% 因此最后乘以 2 得到基带能量估计。
%
% 输入：
%   rx_rf      : 实射频接收信号
%   fs         : 采样率
%   lpf_cutoff : 低通截止频率
%   fir_order  : FIR 滤波器阶数
%
% 输出：
%   energy_est : 估计的 |u(t)|^2

    rx_rf = rx_rf(:);

    rx_square = rx_rf.^2;

    b_lpf = fir1(fir_order, lpf_cutoff/(fs/2), 'low');

    rx_energy_lpf = filtfilt(b_lpf, 1, rx_square);

    energy_est = 2 * rx_energy_lpf;

    energy_est = max(real(energy_est), 0);

    energy_est = energy_est(:);
end


function E_symbol = extract_symbol_energy_from_waveform(energy_waveform, sym_center_idx, win_len)
%EXTRACT_SYMBOL_ENERGY_FROM_WAVEFORM 从能量波形中提取每个符号的平均能量
%
% 输入：
%   energy_waveform : 平方律检波后的能量波形
%   sym_center_idx  : 每个符号中心索引
%   win_len         : 每个符号窗口长度
%
% 输出：
%   E_symbol        : 1 x N_sym 的符号平均能量

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


function Cp = estimate_pulse_energy_coefficient( ...
    N_sym, rrc, sps, interp, sym_center_idx, win_len)
%ESTIMATE_PULSE_ENERGY_COEFFICIENT 估计符号窗口内的脉冲平均能量系数
%
% 用一个单位符号通过同样的 RRC 成形和插值，计算其在符号窗口内的平均能量。
%
% 输出：
%   Cp : 符号窗口内脉冲平均能量系数

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


function symbols_est = symbol_energy_ml_detector( ...
    E_simo, h_ris, ref_symbols_ris, constellation, Cp)
%SYMBOL_ENERGY_ML_DETECTOR 基于符号级能量的 SIMO ML 检测器
%
% 简化预测模型：
%   E_{m,k} ≈ Cp * | h_m * c + b_{m,k} |^2
%
% 对每个符号 k，遍历星座点 c，选择预测能量与观测能量最接近的星座点。
%
% 输入：
%   E_simo          : M_ris x N_sym 符号级能量观测
%   h_ris           : M_ris x 1 阵列流形
%   ref_symbols_ris : N_sym x M_ris 参考符号矩阵
%   constellation   : 星座点列向量
%   Cp              : 脉冲能量系数
%
% 输出：
%   symbols_est     : N_sym x 1 估计符号

    h_ris = h_ris(:);
    constellation = constellation(:);

    [M_ris, N_sym] = size(E_simo);

    if length(h_ris) ~= M_ris
        error('h_ris 的长度必须等于 E_simo 的行数。');
    end

    if size(ref_symbols_ris,1) ~= N_sym || size(ref_symbols_ris,2) ~= M_ris
        error('ref_symbols_ris 的尺寸必须为 N_sym x M_ris。');
    end

    symbols_est = zeros(N_sym, 1);

    for k = 1:N_sym

        metric = zeros(length(constellation), 1);

        for ci = 1:length(constellation)

            c = constellation(ci);

            E_pred = zeros(M_ris, 1);

            for m = 1:M_ris
                b_mk = ref_symbols_ris(k,m);

                E_pred(m) = Cp * abs(h_ris(m) * c + b_mk).^2;
            end

            metric(ci) = sum((E_simo(:,k) - E_pred).^2);
        end

        [~, id_min] = min(metric);

        symbols_est(k) = constellation(id_min);
    end
end