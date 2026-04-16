%% FMCW 符号能量检测仿真 这个脚本使用了chirp信号作为符号，成功，但是不可以作为后续的基础，只能另想idea
clear; clc; close all;

%% 1. 基本参数
N_sym = 8;              % 符号数
Rs = 1e6;               % 符号速率 1 MHz
Ts = 1/Rs;              % 符号周期 1 μs
sps = 20;               % 每符号采样点数
fs = sps * Rs;          % 采样率 20 MHz
fc = 3.5e9;             % 载波频率 3.5 GHz

% FMCW 参数
B = 5e6;                % FMCW 带宽 5 MHz
T_chirp = Ts;           % 一个符号 = 一个 chirp
k = B / T_chirp;        % 调频斜率

%% 2. 生成 4 种不同的 FMCW 符号（斜率不同）
% 用 4 种斜率代表 QPSK 的 4 个符号
slopes = [-k, -k/3, k/3, k];   % 四种调频斜率

% 随机选择 8 个符号
symbol_idx = randi([1 4], 1, N_sym);

%% 3. 生成 FMCW 基带信号
t_sym = (0:sps-1)/fs;   % 一个符号的时间轴
tx_bb = [];

for i = 1:N_sym
    slope = slopes(symbol_idx(i));
    
    % 线性调频信号（复基带）
    chirp_i = exp(1j*pi*slope*t_sym.^2);
    
    tx_bb = [tx_bb, chirp_i];
end

% 时间轴
t_bb = (0:length(tx_bb)-1)/fs;

%% 4. 上变频到 3.5 GHz
carrier = exp(1j*2*pi*fc*t_bb);
tx_rf = real(tx_bb .* carrier);

%% 5. 能量检测（每个符号周期内）
energy_detected = zeros(1, N_sym);

for i = 1:N_sym
    start_idx = (i-1)*sps + 1;
    end_idx   = i*sps;
    
    segment = tx_rf(start_idx:end_idx);
    energy_detected(i) = sum(abs(segment).^2);
end

%% 6. 可视化
figure('Position', [100 100 1400 800]);

% 子图1：基带 FMCW 时域（实部）
subplot(3,2,1);
plot(t_bb*1e6, real(tx_bb));
xlabel('时间 (\mus)');
ylabel('幅度');
title('基带 FMCW 时域波形（实部）');
grid on;

% 子图2：射频时域
subplot(3,2,2);
plot(t_bb*1e6, tx_rf);
xlabel('时间 (\mus)');
ylabel('幅度');
title('射频信号（3.5 GHz 载波）');
grid on;

% 子图3：基带频谱
subplot(3,2,3);
NFFT = 2^nextpow2(length(tx_bb));
f = (-NFFT/2:NFFT/2-1)*(fs/NFFT)/1e6;
X = fftshift(fft(tx_bb, NFFT));
plot(f, 20*log10(abs(X)/max(abs(X))));
xlabel('频率 (MHz)');
ylabel('归一化幅度 (dB)');
title('基带频谱（带宽可控）');
grid on;

% 子图4：射频频谱
subplot(3,2,4);
X_rf = fftshift(fft(tx_rf, NFFT));
plot(f, 20*log10(abs(X_rf)/max(abs(X_rf))));
xlabel('频率 (MHz)');
ylabel('归一化幅度 (dB)');
title('射频频谱（中心 3500 MHz）');
grid on;

% 子图5：能量检测结果
subplot(3,2,5);
bar(1:N_sym, energy_detected);
xlabel('符号索引');
ylabel('能量');
title('每个符号的能量检测');
grid on;

% 子图6：符号类型
subplot(3,2,6);
stem(1:N_sym, symbol_idx, 'filled');
set(gca, 'YTick', 1:4);
ylabel('符号类型 (1~4)');
xlabel('符号索引');
title('发送的 FMCW 符号类型');
grid on;

%% 7. 输出结果
fprintf('========== FMCW 能量检测仿真 ==========\n');
fprintf('符号数: %d\n', N_sym);
fprintf('符号速率: %.2f MHz\n', Rs/1e6);
fprintf('FMCW 带宽: %.2f MHz\n', B/1e6);
fprintf('载波频率: %.2f GHz\n', fc/1e9);
fprintf('\n每个符号的能量：\n');
disp(energy_detected);
fprintf('能量标准差: %.4f\n', std(energy_detected));
fprintf('能量变化率: %.2f %%\n', std(energy_detected)/mean(energy_detected)*100);