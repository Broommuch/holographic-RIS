%% 标准 QPSK 能量检测（匹配滤波法） 这个脚本用匹配滤波实现了符号的能量检测，对应的也很准，但是不能作为后续的参考
clear; clc; close all;

%% 1. 参数设置
N_sym = 8;              % 8个符号
Rs = 1e6;               % 符号速率
sps = 20;               % 每符号采样点数
fs = sps * Rs;          % 采样率
rolloff = 0.25;         % LTE 标准滚降系数

%% 2. 生成 QPSK 符号
bits = randi([0 1], 1, N_sym*2);
symbols = qammod(bits, 4, 'gray', 'UnitAveragePower', true);

%% 3. 成形滤波（发射端）
rrc = rcosdesign(rolloff, 6, sps, 'sqrt');
tx_bb = upfirdn(symbols, rrc, sps, 1);

% 滤波器延迟
delay = (length(rrc)-1)/2;

%% 4. 接收端：匹配滤波（关键步骤！）
rx_mf = filter(rrc, 1, tx_bb);

% 最佳采样点（符号中心）
sample_points = delay + 1 : sps : length(rx_mf);

% 提取符号
rx_sym = rx_mf(sample_points);

% 计算每个符号的能量
symbol_energy = abs(rx_sym).^2;
symbol_energy = symbol_energy(4:19);

%% 5. 验证：能量是否与发送符号一致？
tx_energy = abs(symbols).^2;
energy_error = mean(abs(symbol_energy - tx_energy)./tx_energy)*100;

fprintf('========== 能量检测结果 ==========\n');
fprintf('符号索引\t发送能量\t接收能量\t误差(%)\n');
for i = 1:N_sym
    fprintf('%d\t\t%.4f\t\t%.4f\t\t%.2f\n', ...
            i, tx_energy(i), symbol_energy(i), ...
            abs(symbol_energy(i)-tx_energy(i))/tx_energy(i)*100);
end
fprintf('\n平均能量误差: %.4f %%\n', energy_error);

%% 6. 时域可视化（关键！）
t = (0:length(tx_bb)-1)/fs * 1e6; % 微秒

figure('Position', [100 100 1200 400]);

% 子图1：发射波形
subplot(1,2,1);
plot(t, real(tx_bb));
xlabel('时间 (\mus)');
ylabel('幅度');
title('发射端成形滤波后（符号边界消失）');
grid on;

% 标记采样点
hold on;
sample_times = (sample_points-1)/fs * 1e6;
stem(sample_times, ones(size(sample_times))*max(real(tx_bb)), ...
     'r.', 'MarkerSize', 10);
legend('波形', '最佳采样点');

% 子图2：能量检测
subplot(1,2,2);
bar(1:N_sym*2, symbol_energy);
xlabel('符号索引');
ylabel('能量');
title('每个符号的能量（匹配滤波后）');
grid on;

%% 7. 如果你想“强行”时域能量检测（非相干）
% 注意：这不如匹配滤波准确，但可用于盲检测
window_len = sps;  % 一个符号长度的窗口
energy_naive = zeros(1, N_sym);

for i = 1:N_sym
    start_idx = (i-1)*sps + 1;
    end_idx = i*sps;
    segment = tx_bb(start_idx:end_idx);
    energy_naive(i) = sum(abs(segment).^2);
end

fprintf('\n========== 非相干能量检测（不准确） ==========\n');
disp(energy_naive);