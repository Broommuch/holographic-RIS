clc; clear; close all;
%% MATLAB脚本：500MHz采样率下对20MHz和1MHz信号采样1000个点

% 参数设置
fs = 500e6;           % 采样率：500 MHz
N = 1000;             % 采样点数
t = (0:N-1)/fs;       % 时间向量

% 信号1：20MHz信号（假设为正弦波）
f1 = 20e6;           % 信号频率：200 MHz
signal1 = 10*sin(2*pi*f1*t);

% 信号2：1MHz正弦信号
f2 = 1e6;             % 信号频率：1 MHz
signal2 = sin(2*pi*f2*t);

signal3 = signal1 .* signal2;

% 创建图形窗口
figure('Position', [100, 100, 1200, 800]);

% 绘制第一个信号
subplot(2,2,1);
plot(t*1e9, signal1, 'b-', 'LineWidth', 1.5);
xlabel('时间 (ns)');
ylabel('幅度');
title('200MHz信号时域波形');
grid on;
xlim([0, 1000]);  % 显示前50ns

subplot(2,2,2);
plot(t*1e9, signal2, 'r-', 'LineWidth', 1.5);
xlabel('时间 (ns)');
ylabel('幅度');
title('1MHz正弦信号时域波形');
grid on;
xlim([0, 1000]);  % 显示前500ns

% 频域分析
subplot(2,2,3);
plot(t*1e9, signal3, 'r-', 'LineWidth', 1.5);
grid on;
xlim([0, 1000]);  % 显示前500ns

% subplot(2,2,4);
% signal2_fft = fftshift(fft(signal2));
% plot(f/1e6, abs(signal2_fft), 'r-', 'LineWidth', 1.5);
% xlabel('频率 (MHz)');
% ylabel('幅度');
% title('1MHz信号频谱');
% grid on;
% xlim([-10, 10]);

% 添加总标题
sgtitle('500MHz采样率下信号采样分析', 'FontSize', 14, 'FontWeight', 'bold');

% 显示采样信息
fprintf('采样参数：\n');
fprintf('采样率: %.2f MHz\n', fs/1e6);
fprintf('采样点数: %d 点\n', N);
fprintf('时间分辨率: %.2f ns\n', 1/fs*1e9);
fprintf('总采样时间: %.2f us\n', N/fs*1e6);
fprintf('\n信号参数：\n');
fprintf('信号1频率: %.2f MHz\n', f1/1e6);
fprintf('信号2频率: %.2f MHz\n', f2/1e6);