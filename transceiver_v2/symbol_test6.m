%% 该脚本实现了时域波形的完整展示，可作为后续符号提取的参考
clc; clear; close all;

%% ================= 参数设置 =================
N_sym = 50;          % 符号数（少一点方便观察）
sps   = 8;           % 每符号采样点（基带）
rolloff = 0.25;
span = 6;

Rs = 1e6;            % 符号率 1 MHz
fs_bb = Rs * sps;    % 基带采样率

interp = 50;         % 额外插值倍数（关键：提高RF分辨率）
fs = fs_bb * interp; % 最终采样率（很高）

fc = 20e6;           % 射频载波（20 MHz）

rng(64);
%% ================= QPSK调制 =================
bits = randi([0 1], 1, 2*N_sym);
bits_reshape = reshape(bits, 2, []).';
symbols_idx = bi2de(bits_reshape, 'left-msb');

symbols = qammod(symbols_idx, 4, 'gray', ...
    'UnitAveragePower', true);

symbols = symbols(:);   % 强制列向量（避免维度坑）

%% ================= RRC成形 =================
rrc = rcosdesign(rolloff, span, sps, 'sqrt');

tx_bb = upfirdn(symbols, rrc, sps, 1) + 1;  % 基带信号

%% ================= 高采样率插值 =================
tx_bb_hi = resample(tx_bb, interp, 1);  % 提高采样率

%% ================= 时间轴 =================
t = (0:length(tx_bb_hi)-1)' / fs;

%% ================= 上变频 =================
carrier = exp(1j*2*pi*fc*t);
tx_rf = real(tx_bb_hi .* carrier);

snr = 50;
tx_rf_noise = awgn(tx_rf,snr);

%% ================= 包络（用于对比） =================
envelope = abs(tx_bb_hi);

%% ================= 画图 =================
figure('Position',[100,100,1200,800]);

% ---- 基带 I/Q ----
subplot(3,1,1);
plot(t*1e6, real(tx_bb_hi), 'b'); hold on;
plot(t*1e6, imag(tx_bb_hi), 'r');
title('成形滤波后的基带信号 (I/Q)');
xlabel('时间 (us)');
ylabel('幅度');
legend('I','Q');
grid on;

% ---- 包络 ----
subplot(3,1,2);
plot(t*1e6, envelope, 'k', 'LineWidth', 1.2);
title('基带包络 |s(t)|');
xlabel('时间 (us)');
ylabel('幅度');
grid on;

% ---- 射频信号 ----
subplot(3,1,3);
plot(t(1:end)*1e6, tx_rf_noise(1:length(t))); % 放大局部
title('射频信号（可见载波 + 包络）');
xlabel('时间 (us)');
ylabel('幅度');
grid on;

sgtitle('QPSK + RRC 成形 + 射频调制（高采样率可视化）');

%% ================= 信息输出 =================
fprintf('======= 参数信息 =======\n');
fprintf('符号数: %d\n', N_sym);
fprintf('基带采样率: %.2f MHz\n', fs_bb/1e6);
fprintf('最终采样率: %.2f MHz\n', fs/1e6);
fprintf('载波频率: %.2f MHz\n', fc/1e6);
fprintf('总采样点数: %d\n', length(tx_bb_hi));

%% ================= 接收端（非相干能量检测） =================

% ---- 包络提取 ----
rx_env = abs(hilbert(tx_rf));

%% ===== 去除滤波器延迟 =====
delay = span * sps / 2 * interp;

N = length(rx_env);
rx_env_valid = rx_env(delay+1 : N - delay);

%% ===== 符号能量提取 =====
samples_per_symbol = sps * interp;

N_sym_rx = floor(length(rx_env_valid)/samples_per_symbol);

energy = zeros(N_sym_rx,1);

for k = 1:N_sym_rx
    idx = (k-1)*samples_per_symbol + (1:samples_per_symbol);
    energy(k) = sum(rx_env_valid(idx).^2);
end

%% ===== 可视化 =====
figure;
stem(energy, 'filled');
title(['检测到符号数 = ', num2str(N_sym_rx)]);
ylabel('energy of each symbol');
grid on;