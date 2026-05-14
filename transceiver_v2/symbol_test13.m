% 这个脚本尝试改变数据模型为simo，重新梳理代码逻辑 
% 改成了simo模型，但是是在波形域进行的检测，不是在能量积分域，再看看能不能在能量积分域检测

clc; clear; close all;

%% ================= 参数设置 =================
N_sym = 20;          % 符号数（少一点方便观察）
N_sym_ref = N_sym;   % 参考信号符号数等于未知信号符号数
sps   = 8;           % 每符号采样点（基带）
rolloff = 0.25;
span = 6;

Rs = 1e6;            % 符号率 1 MHz
fs_bb = Rs * sps;    % 基带采样率

interp = 50;         % 额外插值倍数（关键：提高RF分辨率）
fs = fs_bb * interp; % 最终采样率（很高）

fc = 20e6;           % 射频载波（20 MHz）
fc_ref = 20e6;       % 参考信号频率

rng(64);
%% ================= QPSK调制 =================
bits = randi([0 1], 1, 2*N_sym);
bits_reshape = reshape(bits, 2, []).';
symbols_idx = bi2de(bits_reshape, 'left-msb');

symbols = qammod(symbols_idx, 4, 'gray', ...
    'UnitAveragePower', true);

% symbols = exp(1i*2*pi*rand(N_sym, 1));  % True symbol vector (complex unit modulus) 尝试先更改发送符号
symbols = symbols(:);   % 强制列向量（避免维度坑）
% symbols_ref = ones(N_sym_ref,1)*symbols(1); % 参考信号符号不变

%% ================= 参考符号构造方法1：构造随机相位参考符号序列 =================

ref_amp = 1.5;

rng(2026);
ref_phase_seq = 2*pi*rand(N_sym_ref, 1);

symbols_ref = ref_amp * exp(1j * ref_phase_seq);
ref_symbol = symbols_ref;

%% =================  参考符号构造方法2：构造确定性变化参考序列 =================

ref_amp = 1.5;

k_ref = (0:N_sym_ref-1).';

ref_phase_seq = mod(2*pi*0.137*k_ref.^2 + pi/7*k_ref, 2*pi);

symbols_ref = ref_amp * exp(1j * ref_phase_seq);

%% ================= RRC成形 =================
rrc = rcosdesign(rolloff, span, sps, 'sqrt');

tx_bb = upfirdn(symbols, rrc, sps, 1) ;  % 基带信号
tx_bb_ref = upfirdn(symbols_ref, rrc, sps, 1) ;  % 参考信号的基带信号

%% ================= 高采样率插值 =================
tx_bb_hi = resample(tx_bb, interp, 1);  % 提高采样率
tx_bb_hi_ref = resample(tx_bb_ref, interp, 1);  % 提高采样率

%% ================= 时间轴 =================
t = (0:length(tx_bb_hi)-1)' / fs;

%% ================= 上变频 =================
carrier = exp(1j*2*pi*fc*t);
tx_rf = real(tx_bb_hi .* carrier);

snr = 20;
tx_rf_noise = awgn(tx_rf,snr);

%% ================= 叠加参考信号===============
% phi_ref = -1*pi/4;
% carrier_ref = exp(1j*2*pi*fc*t + phi_ref);
tx_rf_ref = real(tx_bb_hi_ref .* carrier);
ref_rf = real(tx_rf_ref);
mixed = tx_rf + ref_rf;

mixed_rf_noise = awgn(mixed,snr);


%% ================= 包络（用于对比） =================
envelope = abs(tx_bb_hi);
envelope_mixed = abs(mixed_rf_noise);

%% ================= 画图 =================
figure('Position',[100,100,1200,800]);

fig_num = 5;

% ---- 基带 I/Q ----
subplot(fig_num,1,1);
plot(t*1e6, real(tx_bb_hi), 'b'); hold on;
plot(t*1e6, imag(tx_bb_hi), 'r');
title('成形滤波后的基带信号 (I/Q)');
xlabel('时间 (us)');
ylabel('幅度');
legend('I','Q');
grid on;

% ---- 包络 ----
subplot(fig_num,1,2);
plot(t*1e6, envelope, 'k', 'LineWidth', 1.2);
title('基带包络 |s(t)|');
xlabel('时间 (us)');
ylabel('幅度');
grid on;

% ---- 射频信号 ----
subplot(fig_num,1,3);
plot(t(1:end)*1e6, tx_rf_noise(1:length(t))); % 放大局部
title('射频信号（可见载波 + 包络）');
xlabel('时间 (us)');
ylabel('幅度');
grid on;

% ---- 参考信号 ----
subplot(fig_num,1,4);
plot(t(1:end)*1e6, ref_rf(1:length(t))); % 放大局部
title('参考信号');
xlabel('时间 (us)');
ylabel('幅度');
grid on;

% ---- 叠加后的射频信号 ----
subplot(fig_num,1,5);
plot(t(1:end)*1e6, mixed_rf_noise(1:length(t))); % 放大局部
title('叠加后的射频信号');
xlabel('时间 (us)');
ylabel('幅度');
grid on;

sgtitle('QPSK + RRC 成形 + 射频调制（高采样率可视化）');

%% ================= 平方律检波 + 低通滤波 =================

% 选择要检测的实通带信号
% 如果想看无噪声情况，用 mixed
% 如果想看有噪声情况，用 mixed_rf_noise
rx_rf = mixed_rf_noise(:);

% 1. 平方律检波
rx_square = rx_rf.^2;

% 2. 设计低通滤波器
% 基带信号带宽大约为 Rs*(1+rolloff)，平方后包络带宽会扩展一些
% 这里取 5 MHz，远小于 2fc = 40 MHz，可以滤掉二倍载频项
lpf_cutoff = 5e6;          % 低通截止频率
fir_order = 800;           % FIR阶数，可根据平滑程度调整

b_lpf = fir1(fir_order, lpf_cutoff/(fs/2), 'low');

% 3. 零相位低通滤波，避免群时延
rx_energy_lpf = filtfilt(b_lpf, 1, rx_square);

% 4. 根据 s_RF(t)=real{z_bb(t) exp(j2πfct)} 的定义，需要乘以2
rx_energy_est = 2 * rx_energy_lpf;

%% ================= 理论复基带能量用于对比 =================

% 因为你的两个实通带信号同频同相叠加：
% real(tx_bb_hi*carrier) + real(tx_bb_hi_ref*carrier)
% = real((tx_bb_hi + tx_bb_hi_ref)*carrier)
z_bb_hi = tx_bb_hi(:) + tx_bb_hi_ref(:);

energy_theory = abs(z_bb_hi).^2;

%% ================= 去除滤波边缘过渡区 =================

trim = span*sps*interp + fir_order;

valid_idx = (trim+1):(length(rx_energy_est)-trim);

t_valid = t(valid_idx);
energy_est_valid = rx_energy_est(valid_idx);
energy_theory_valid = energy_theory(valid_idx);

%% ================= 结果可视化 =================

figure;
plot(t_valid*1e6, energy_theory_valid, 'LineWidth', 1.5); hold on;
plot(t_valid*1e6, energy_est_valid, '--', 'LineWidth', 1.2);
grid on;
xlabel('Time / \mus');
ylabel('Energy');
legend('Theoretical |z_{bb}(t)|^2', 'Detected energy after square-law + LPF');
title('Square-law Detection Result');

%% ================= 局部放大观察 =================

figure;
idx_show = valid_idx(1:min(4000, length(valid_idx)));

plot(t(idx_show)*1e6, energy_theory(idx_show), 'LineWidth', 1.5); hold on;
plot(t(idx_show)*1e6, rx_energy_est(idx_show), '--', 'LineWidth', 1.2);
grid on;
xlabel('Time / \mus');
ylabel('Energy');
legend('Theoretical |z_{bb}(t)|^2', 'Detected energy');
title('Zoomed View of Detected Envelope Energy');

%% ================= 误差评估 =================

mse_energy = mean(abs(energy_est_valid - energy_theory_valid).^2);
nmse_energy = mse_energy / mean(abs(energy_theory_valid).^2);

fprintf('Energy detection MSE  = %.4e\n', mse_energy);
fprintf('Energy detection NMSE = %.4e\n', nmse_energy);

%% ================= 符号级能量切分 =================

% 每个符号对应的高采样率采样点数
sym_samp_hi = sps * interp;

% RRC滤波器群延迟
% rrc长度为 span*sps + 1
% 群延迟为 span*sps/2 个基带采样点
rrc_delay_bb = span * sps / 2;

% 转换到高采样率后的群延迟
rrc_delay_hi = rrc_delay_bb * interp;

% 每个符号中心位置
% MATLAB索引从1开始
sym_center_idx = round(rrc_delay_hi + 1 + (0:N_sym-1).' * sym_samp_hi);

% 预分配
symbol_energy_center = zeros(N_sym, 1);   % 符号中心采样能量
symbol_energy_avg    = zeros(N_sym, 1);   % 符号窗口平均能量
symbol_energy_int    = zeros(N_sym, 1);   % 符号窗口积分能量

% 窗口长度：一个符号周期
win_len = sym_samp_hi;

for k = 1:N_sym

    % 当前符号中心
    cidx = sym_center_idx(k);

    % 方法1：符号中心采样
    if cidx >= 1 && cidx <= length(rx_energy_est)
        symbol_energy_center(k) = rx_energy_est(cidx);
    else
        symbol_energy_center(k) = NaN;
    end

    % 方法2：以符号中心为中心，取一个符号周期窗口
    left_idx  = round(cidx - win_len/2);
    right_idx = round(cidx + win_len/2 - 1);

    % 防止越界
    left_idx  = max(left_idx, 1);
    right_idx = min(right_idx, length(rx_energy_est));

    % 当前符号窗口内的能量波形
    energy_window = rx_energy_est(left_idx:right_idx);

    % 平均能量
    symbol_energy_avg(k) = mean(energy_window);

    % 积分能量
    % 离散积分需要乘以采样间隔 1/fs
    symbol_energy_int(k) = sum(energy_window) / fs;
end

%% ================= 理论符号级能量对比 =================

% 理论连续能量
energy_theory = abs(z_bb_hi).^2;

symbol_energy_theory_center = zeros(N_sym, 1);
symbol_energy_theory_avg    = zeros(N_sym, 1);
symbol_energy_theory_int    = zeros(N_sym, 1);

for k = 1:N_sym

    cidx = sym_center_idx(k);

    if cidx >= 1 && cidx <= length(energy_theory)
        symbol_energy_theory_center(k) = energy_theory(cidx);
    else
        symbol_energy_theory_center(k) = NaN;
    end

    left_idx  = round(cidx - win_len/2);
    right_idx = round(cidx + win_len/2 - 1);

    left_idx  = max(left_idx, 1);
    right_idx = min(right_idx, length(energy_theory));

    energy_window_theory = energy_theory(left_idx:right_idx);

    symbol_energy_theory_avg(k) = mean(energy_window_theory);
    symbol_energy_theory_int(k) = sum(energy_window_theory) / fs;
end

%% ================= 可视化：符号级能量序列 =================

figure;
stem(1:N_sym, symbol_energy_theory_avg, 'LineWidth', 1.5); hold on;
stem(1:N_sym, symbol_energy_avg, '--', 'LineWidth', 1.2);
grid on;
xlabel('Symbol index');
ylabel('Average energy');
legend('Theoretical symbol energy', 'Detected symbol energy');
title('Symbol-level Energy after Square-law Detection');

%% ================= 可视化：连续能量与符号窗口中心 =================

figure;
plot(t*1e6, energy_theory, 'LineWidth', 1.2); hold on;
plot(t*1e6, rx_energy_est, '--', 'LineWidth', 1.0);

valid_center_idx = sym_center_idx;
valid_center_idx = valid_center_idx(valid_center_idx >= 1 & valid_center_idx <= length(t));

scatter(t(valid_center_idx)*1e6, rx_energy_est(valid_center_idx), 35, 'filled');

grid on;
xlabel('Time / \mus');
ylabel('Energy');
legend('Theoretical |z_{bb}(t)|^2', 'Detected energy', 'Symbol centers');
title('Continuous Energy Waveform and Symbol Centers');

%% ================= 输出最终的50个符号能量 =================

% 推荐后续算法使用这个：
E_symbol = symbol_energy_avg;

disp('50个符号对应的混合参考信号后的平均能量为：');
disp(E_symbol.');

%% ================= 误差评估 =================

symbol_nmse = mean(abs(symbol_energy_avg - symbol_energy_theory_avg).^2) ...
              / mean(abs(symbol_energy_theory_avg).^2);

fprintf('Symbol-level energy NMSE = %.4e\n', symbol_nmse);

%% ================= RRC波形域GS检测：构造 P 矩阵 =================

% ================= 增加观测点：每个符号周期内取更多采样点 =================

sym_samp_hi = sps * interp;
rrc_delay_hi = span * sps / 2 * interp;

% 每个符号内取 M_obs_per_sym 个观测点
% 可以先试 21，再试 41
M_obs_per_sym = 11;

% 在每个符号周期内均匀取点，避开边界附近
obs_offsets = round(linspace(-0.48, 0.48, M_obs_per_sym) * sym_samp_hi);
obs_offsets = unique(obs_offsets);

obs_idx = [];

for k = 1:N_sym
    center_k = round(rrc_delay_hi + 1 + (k-1) * sym_samp_hi);

    for q = 1:length(obs_offsets)
        idx = center_k + obs_offsets(q);

        if idx >= 1 && idx <= length(tx_bb_hi)
            obs_idx = [obs_idx; idx];
        end
    end
end

obs_idx = unique(obs_idx);
N_obs = length(obs_idx);

fprintf('Dense observation mode: N_obs = %d, N_sym = %d, ratio = %.2f\n', ...
        N_obs, N_sym, N_obs/N_sym);

%% ================= 构造 P: symbols -> shaped waveform samples =================

P = zeros(N_obs, N_sym);

for m = 1:N_sym

    s_basis = zeros(N_sym, 1);
    s_basis(m) = 1;

    bb_basis = upfirdn(s_basis, rrc, sps, 1);
    bb_basis_hi = resample(bb_basis, interp, 1);

    for n = 1:N_obs
        idx = obs_idx(n);

        if idx >= 1 && idx <= length(bb_basis_hi)
            P(n, m) = bb_basis_hi(idx);
        else
            P(n, m) = 0;
        end
    end
end


%% ================= SIMO-RIS 阵列参数 =================

RIS_row = 8;
RIS_col = 8;
M_ris = RIS_row * RIS_col;

lambda = 1;          % 归一化波长
d = 0.5 * lambda;    % 阵元间距
k0 = 2*pi/lambda;

% 单用户到达角，按你的定义调整
theta_u = 20 * pi/180;   % 方位/俯仰需要和你的坐标系保持一致
phi_u   = 10 * pi/180;

% 构造 RIS 单元坐标，这里假设 RIS 位于 y-z 平面，x 为法向
% row 对应 z 方向，col 对应 y 方向
[y_idx, z_idx] = meshgrid(0:RIS_col-1, 0:RIS_row-1);

y_pos = (y_idx(:) - (RIS_col-1)/2) * d;
z_pos = (z_idx(:) - (RIS_row-1)/2) * d;

% 入射方向在 y-z 平面的方向余弦
% 这里沿用你之前常用的形式：
% ky = k sin(theta) cos(phi)
% kz = k sin(theta) sin(phi)
ky = k0 * sin(theta_u) * cos(phi_u);
kz = k0 * sin(theta_u) * sin(phi_u);

% 单用户到每个 RIS 单元的阵列流形
a_ris = exp(1j * (ky * y_pos + kz * z_pos));   % M_ris x 1

% 可选：加入每个单元的幅度增益
% amp_ris = ones(M_ris,1);
% a_ris = amp_ris .* a_ris;

%% ================= 每个 RIS 单元的独立参考符号 =================

ref_amp = 1.5;

% 基础参考序列，可以仍然用你之前的确定性变化序列
k_ref = (0:N_sym-1).';
base_ref_phase = mod(2*pi*0.137*k_ref.^2 + pi/7*k_ref, 2*pi);

% 每个 RIS 单元额外加一个独立参考相位
rng(2028);
ris_ref_phase = 2*pi*rand(M_ris, 1);

symbols_ref_ris = zeros(N_sym, M_ris);

for m = 1:M_ris
    symbols_ref_ris(:,m) = ref_amp * exp(1j * (base_ref_phase + ris_ref_phase(m)));
end

%% ================= 构造 SIMO 等效矩阵 =================
% 单通道模型：
% z_m = | a_m * P * s + P * b_m |
%
% 堆叠模型：
% z_all = | P_simo * s + b_simo |

P_simo = zeros(M_ris * N_obs, N_sym);
b_simo = zeros(M_ris * N_obs, 1);

for m = 1:M_ris
    rows = (m-1)*N_obs + (1:N_obs);

    P_simo(rows, :) = a_ris(m) * P;
    b_simo(rows) = P * symbols_ref_ris(:,m);
end

%% ================= 生成 SIMO 理论观测幅值 =================

z_simo = abs(P_simo * symbols(:) + b_simo);

fprintf('SIMO observation size: %d x 1\n', length(z_simo));
fprintf('Unknown symbol size  : %d x 1\n', N_sym);
fprintf('Observation ratio    : %.2f\n', length(z_simo)/N_sym);

%% ================= SIMO 前向模型一致性检查 =================

z_model_simo = abs(P_simo * symbols(:) + b_simo);

simo_mag_error = norm(z_simo - z_model_simo) / norm(z_model_simo);

fprintf('SIMO magnitude model relative error = %.4e\n', simo_mag_error);

figure;
plot(z_model_simo, 'LineWidth', 1.2); hold on;
plot(z_simo, '--', 'LineWidth', 1.0);
grid on;
xlabel('Stacked observation index');
ylabel('Magnitude');
legend('Theoretical SIMO |A s + b|', 'Observed SIMO magnitude');
title('SIMO Magnitude Observation Check');

% ---------- 调用恢复算法 ----------
t0 = 1000;

% 推荐使用统一接口：z = abs(A*s + b)
s_est = biased_gn_algorithm(z_simo, P_simo, b_simo, t0);

% % 如果还想调用原GS接口：
% A_gs_simo = P_simo.';
% s_est = biased_gs_algorithm(z_simo, A_gs_simo, b_simo, t0);

%% ================= QPSK硬判决 =================

qpsk_const = qammod((0:3).', 4, 'gray', 'UnitAveragePower', true);

s_detect = zeros(N_sym, 1);

for k = 1:N_sym
    [~, idx_min] = min(abs(s_est(k) - qpsk_const));
    s_detect(k) = qpsk_const(idx_min);
end

symbol_error = sum(s_detect ~= symbols(:));
SER = symbol_error / N_sym;

fprintf('Waveform-domain GS symbol errors = %d / %d\n', symbol_error, N_sym);
fprintf('Waveform-domain GS SER = %.4f\n', SER);

%% ================= 可视化恢复结果 =================

figure;
plot(real(symbols), imag(symbols), 'o', 'LineWidth', 1.5); hold on;
plot(real(s_est), imag(s_est), 'x', 'LineWidth', 1.5);
plot(real(s_detect), imag(s_detect), 's', 'LineWidth', 1.2);
grid on;
axis equal;
xlabel('In-phase');
ylabel('Quadrature');
legend('True QPSK symbols', 'GS estimated symbols', 'Hard-decided symbols');
title('Waveform-domain GS Recovery with RRC Shaping');