%% ================= 简陋射频接收机仿真 =================
clear; clc; close all;

%% 1. 参数设置
N_sym = 8;
Rs = 1e6;          % 符号速率 1 MHz
sps = 20;          % 每符号采样点数
fs = sps * Rs;     % 采样率 20 MHz
fc = 3.5e9;        % 载波频率 (注意：这里用基带等效，实际射频需更高采样率)
rolloff = 0.25;

%% 2. 发射端（标准 QPSK + RRC）
bits = randi([0 1], 1, N_sym*2);
symbols = qammod(bits, 4, 'gray', 'UnitAveragePower', true);

rrc = rcosdesign(rolloff, 6, sps, 'sqrt');
tx_bb = upfirdn(symbols, rrc, sps, 1);

% figure;
% plot(tx_bb);

% 上变频（模拟射频）
t = (0:length(tx_bb)-1)/fs;
tx_rf = real(tx_bb .* exp(1j*2*pi*fc*t));
figure;
plot(tx_rf);
%% 3. 简陋接收机处理（无下变频，无匹配滤波）
% 3.1 平方律检波（模拟二极管/包络检波）
P_rf = abs(tx_rf).^2;  % 瞬时功率

% 3.2 滑动窗口积分（关键步骤）
% 积分窗口长度 = 一个符号的采样点数
window_len = sps; 
integral_energy = filter(ones(1, window_len), 1, P_rf);

% 3.3 在符号中心取能量值（对齐符号）
delay = floor(window_len/2);
energy_detected = integral_energy(delay + 1 : sps : end);

% 去掉可能多出来的点
energy_detected = energy_detected(1:N_sym);

%% 4. 结果分析
fprintf('========== 简陋接收机能量检测 ==========\n');
fprintf('符号索引\t检测能量\t实际符号\n');
for i = 1:N_sym
    fprintf('%d\t\t%.4f\t\t%d\n', i, energy_detected(i), bits(2*i-1)*2+bits(2*i));
end

%% 5. 绘图
figure('Position', [100 100 1000 600]);

subplot(3,1,1);
plot(t*1e6, tx_rf);
xlabel('时间 (\mus)');
ylabel('幅度');
title('射频信号（含载波）');
grid on;

subplot(3,1,2);
plot(t*1e6, P_rf);
xlabel('时间 (\mus)');
ylabel('功率');
title('平方律检波后的功率（无载波）');
grid on;

subplot(3,1,3);
stem(1:N_sym, energy_detected, 'filled');
xlabel('符号索引');
ylabel('积分能量');
title('滑动窗口积分后的符号能量');
grid on;