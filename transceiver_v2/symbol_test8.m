%% 这个脚本准备修复gs检测阶段出现的问题

clc; clear; close all;

%% ================= 参数设置 =================
N_sym = 50;          % 符号数（少一点方便观察）
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

% 更改参考信号为已知固定符号
ref_phase = pi/7;                 % 不要取 0, pi/4, pi/2 这类对称角
ref_amp = 1.2;                    % 参考幅度可以略大于信号幅度
ref_symbol = ref_amp * exp(1j*ref_phase);

qpsk_const = qammod((0:3).', 4, 'gray', 'UnitAveragePower', true);

E_qpsk_ref = abs(qpsk_const + ref_symbol).^2;

disp('四个QPSK点叠加参考后的理论能量：');
disp(E_qpsk_ref.');

% 验证四个符号能量是否可区分
figure;
stem(0:3, E_qpsk_ref, 'LineWidth', 1.5);
grid on;
xlabel('QPSK symbol index');
ylabel('|s + b|^2');
title('Energy Separability under Current Reference Symbol');

symbols_ref = ones(N_sym_ref,1) * ref_symbol;

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

snr = 50;
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

%% ================= 用真实参考信号构造 GS 检测模型 =================

% 观测使用符号中心能量
E_obs = symbol_energy_center(:);
E_obs = max(E_obs, 0);
z_obs = sqrt(E_obs);

K_gs = N_sym;       % 待恢复的QPSK符号数量
N_obs = N_sym;      % 目前每个符号中心一个观测

%% ================= 构造 A_eff_H，使得 A_eff_H * symbols ≈ shaped_unknown(center) =================

A_eff_H = zeros(N_obs, K_gs);

for m = 1:K_gs

    s_basis = zeros(K_gs, 1);
    s_basis(m) = 1;

    % 只构造未知符号经过RRC成形后的贡献，不加 +1
    bb_basis = upfirdn(s_basis, rrc, sps, 1);
    bb_basis_hi = resample(bb_basis, interp, 1);

    for n = 1:N_obs
        idx = sym_center_idx(n);

        if idx >= 1 && idx <= length(bb_basis_hi)
            A_eff_H(n, m) = bb_basis_hi(idx);
        else
            A_eff_H(n, m) = 0;
        end
    end
end

% biased_gs_algorithm 使用 z = abs(A' * s + b)
A_gs = A_eff_H.';

%% ================= 构造真实一致的 b_eff =================

% 真实叠加基带：
% z_bb_hi = tx_bb_hi + tx_bb_hi_ref
%
% tx_bb_hi = shaped_unknown 
% tx_bb_hi_ref = shaped_reference 
%
% A_gs' * s 表示 shaped_unknown
% 所以 b_eff = tx_bb_hi_ref 

b_eff = zeros(N_obs, 1);

for n = 1:N_obs
    idx = sym_center_idx(n);

    if idx >= 1 && idx <= length(tx_bb_hi_ref)
        b_eff(n) = tx_bb_hi_ref(idx) ;
    else
        b_eff(n) = 1;
    end
end

%% ================= 检查模型是否和真实基带一致 =================

z_model_true = A_gs' * symbols + b_eff;
z_real_true  = zeros(N_obs, 1);

for n = 1:N_obs
    idx = sym_center_idx(n);

    if idx >= 1 && idx <= length(z_bb_hi)
        z_real_true(n) = z_bb_hi(idx);
    else
        z_real_true(n) = NaN;
    end
end

model_error = norm(z_model_true - z_real_true) / norm(z_real_true);

fprintf('Baseband model relative error = %.4e\n', model_error);

figure;
plot(abs(z_real_true), 'o-', 'LineWidth', 1.3); hold on;
plot(abs(z_model_true), 'x--', 'LineWidth', 1.2);
grid on;
xlabel('Symbol index');
ylabel('Magnitude');
legend('Actual |z_{bb}| at centers', 'Model |A^H s + b|');
title('Check of GS Forward Model');

%% ================= 运行 biased GS =================

t0 = 800;

% s_est = biased_gs_algorithm(z_obs, A_gs, b_eff, t0);
s_est = biased_gs_algorithm(z_obs, A_gs, symbols_ref, t0); % 尝试直接调用参考信号估计

%% ================= QPSK硬判决 =================

qpsk_const = qammod((0:3).', 4, 'gray', 'UnitAveragePower', true);

s_detect = zeros(K_gs, 1);
idx_detect = zeros(K_gs, 1);

for k = 1:K_gs
    [~, idx_min] = min(abs(s_est(k) - qpsk_const));
    s_detect(k) = qpsk_const(idx_min);
    idx_detect(k) = idx_min - 1;
end


s_detect = custom_symbol_decision(s_est);
symbol_error = sum(abs(s_detect - symbols) > 0.1);
SER = symbol_error / N_sym;

fprintf('QPSK symbol errors = %d / %d\n', symbol_error, N_sym);
fprintf('SER = %.4f\n', SER);

figure;
plot(real(symbols), imag(symbols), 'o', 'LineWidth', 1.5); hold on;
plot(real(s_est), imag(s_est), 'x', 'LineWidth', 1.5);
plot(real(s_detect), imag(s_detect), 's', 'LineWidth', 1.2);
grid on;
axis equal;
xlabel('In-phase');
ylabel('Quadrature');
legend('True QPSK symbols', 'GS estimated symbols', 'Hard-decided symbols');
title('QPSK Recovery with Correct Reference Field');