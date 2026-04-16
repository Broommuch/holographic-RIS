%%  QPSK 射频能量检测仿真 这个脚本用成形滤波设计了几个符号，但是依旧有很大的问题，不可以直接用，但是让人看到了希望，说不定之后可以做到
clear; clc; close all;

%% 1. 基本参数设置
N_sym = 2;              % 符号数量
Rs = 1e6;               % 符号速率 1 MHz
sps = 20;               % 每符号采样点数
fs = sps * Rs;          % 采样率 20 MHz
fc = 3.5e9;             % 载波频率 3.5 GHz
rolloff = 0.25;         % 滚降系数
span = 6;               % 滤波器跨度

%% 2. 生成 QPSK 符号
% 生成随机比特
bits = randi([0 1], 1, N_sym*2);
bits = [1,1,1,0];
% QPSK 调制（Gray 编码）
symbols = zeros(1, N_sym);
for i = 1:N_sym
    bit_pair = bits(2*i-1:2*i);
    if isequal(bit_pair, [0 0])
        symbols(i) = exp(1j*pi/4);      % 45°
    elseif isequal(bit_pair, [0 1])
        symbols(i) = exp(1j*3*pi/4);    % 135°
    elseif isequal(bit_pair, [1 1])
        symbols(i) = exp(1j*5*pi/4);    % 225°
    else % [1 0]
        symbols(i) = exp(1j*7*pi/4);    % 315°
    end
end

fprintf('生成的 QPSK 符号（复数）:\n');
disp(symbols);

%% 3. 成形滤波（RRC）
rrc = rcosdesign(rolloff, span, sps, 'sqrt');

% 上采样 + 成形滤波
tx_bb = upfirdn(symbols, rrc, sps, 1);
figure;
plot(real(tx_bb));

% 计算群延迟并补偿
delay = (length(rrc)-1)/2;
tx_bb_aligned = tx_bb(delay+1:end-delay);  % 去除延迟

% 归一化能量
tx_bb_norm = tx_bb_aligned / sqrt(mean(abs(tx_bb_aligned).^2));

figure;
plot(real(tx_bb_norm));

%% 4. 上变频到射频
t_bb = (0:length(tx_bb_norm)-1)/fs;
carrier = exp(1j*2*pi*fc*t_bb);
tx_rf = real(tx_bb_norm .* carrier);  % 实部为实信号

%% 5. 射频功率检测（能量检测）
% 创建时间向量
t_rf = (0:length(tx_rf)-1)/fs;

% 方法1：直接对射频信号进行能量检测（滑动窗口积分）
window_len = sps;  % 一个符号周期的采样点数
energy_rf = zeros(1, N_sym);

for i = 1:N_sym
    start_idx = (i-1)*sps + 1;
    end_idx = i*sps;
    
    % 提取一个符号周期的信号
    if i ~= N_sym
        sym_segment = tx_rf(start_idx:end_idx);
    else
        sym_segment = tx_rf(start_idx:end);
    end
    % 计算能量（平方和）
    energy_rf(i) = sum(abs(sym_segment).^2);
end

% 方法2：通过匹配滤波恢复符号能量（参考值）
rx_mf = filter(rrc, 1, tx_bb_norm);
rx_sym = rx_mf(delay+1:sps:end);
energy_ref = abs(rx_sym).^2;

%% 6. 可视化结果
figure('Position', [100 100 1200 800]);

% 子图1：基带时域波形
subplot(3,2,1);
plot(t_bb*1e6, real(tx_bb_norm));
xlabel('时间 (\mus)');
ylabel('幅度');
title('基带 QPSK 时域波形（实部）');
grid on;

% 子图2：射频时域波形
subplot(3,2,2);
plot(t_rf*1e6, tx_rf);
xlabel('时间 (\mus)');
ylabel('幅度');
title('射频信号时域波形（3.5 GHz 载波）');
grid on;

% 子图3：基带频谱
subplot(3,2,3);
NFFT = 2^nextpow2(length(tx_bb_norm));
f = (-NFFT/2:NFFT/2-1)*(fs/NFFT)/1e6;
X = fftshift(fft(tx_bb_norm, NFFT));
plot(f, 20*log10(abs(X)/max(abs(X))));
xlabel('频率 (MHz)');
ylabel('归一化幅度 (dB)');
title('基带信号频谱');
grid on;

% 子图4：射频频谱
subplot(3,2,4);
X_rf = fftshift(fft(tx_rf, NFFT));
plot(f, 20*log10(abs(X_rf)/max(abs(X_rf))));
xlabel('频率 (MHz)');
ylabel('归一化幅度 (dB)');
title('射频信号频谱（中心频率 3500 MHz）');
grid on;

% 子图5：能量检测结果
subplot(3,2,5);
symbol_idx = 1:N_sym;
bar(symbol_idx, energy_rf);
hold on;
plot(symbol_idx, energy_ref, 'ro-', 'LineWidth', 2);
xlabel('符号索引');
ylabel('能量');
title('每个符号的能量检测');
legend('射频直接检测', '匹配滤波参考值', 'Location', 'best');
grid on;

% 子图6：星座图
subplot(3,2,6);
scatter(real(rx_sym), imag(rx_sym), 'filled');
axis square;
grid on;
xlabel('同相分量');
ylabel('正交分量');
title('接收端恢复的星座图');
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);

%% 7. 显示数值结果
fprintf('\n========== 能量检测结果 ==========\n');
fprintf('%-10s %-15s %-15s %-10s\n', '符号', '射频检测能量', '参考能量', '误差(%)');
for i = 1:N_sym
    error_percent = abs(energy_rf(i) - energy_ref(i)) / energy_ref(i) * 100;
    fprintf('%-10d %-15.6f %-15.6f %-10.2f\n', ...
            i, energy_rf(i), energy_ref(i), error_percent);
end

fprintf('\n平均能量误差: %.4f %%\n', mean(abs(energy_rf - energy_ref)./energy_ref)*100);

%% 8. 能量检测性能分析
% 计算信噪比
signal_power = mean(energy_ref);
noise_power = var(tx_rf - mean(tx_rf));  % 简化噪声估计
snr_linear = signal_power / noise_power;
snr_db = 10*log10(snr_linear);

fprintf('\n========== 系统参数 ==========\n');
fprintf('符号速率: %.2f MHz\n', Rs/1e6);
fprintf('采样率: %.2f MHz\n', fs/1e6);
fprintf('每符号采样点数: %d\n', sps);
fprintf('载波频率: %.2f GHz\n', fc/1e9);
fprintf('滚降系数: %.2f\n', rolloff);
fprintf('估计SNR: %.2f dB\n', snr_db);