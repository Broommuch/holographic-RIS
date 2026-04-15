%% Gerchberg-Saxton Algorithm for Holographic Phase Retrieval
clc; clear; close all;

%% 1. 生成示例信号 (模拟通信中的复信号)
rng(1);

% 参数设置
N = 128;           % 阵列大小
d = 0.5;           % 阵元间距（波长单位）
theta = 30;        % 方位角 (度)
phi = 20;          % 俯仰角 (度)

% 转弧度
theta_rad = deg2rad(theta);
phi_rad = deg2rad(phi);

% 阵列索引
[x_idx, y_idx] = meshgrid(0:N-1, 0:N-1);
center_x = (N-1)/2;
center_y = (N-1)/2;

% 方向余弦
u = sin(theta_rad) * cos(phi_rad);
v = sin(theta_rad) * sin(phi_rad);

% 目标相位矩阵
phase_relative = 2*pi * d * ((x_idx - center_x)*u + (y_idx - center_y)*v);

% 幅度
O_amp = ones(N,N) + (2*rand(N)-1)/1000;

% 复信号
O = O_amp .* exp(1i * phase_relative);

% 相位
O_ang = angle(O);


%%
% N = 128; % 信号尺寸
% % 生成二维 QAM 信号图像
[x, y] = meshgrid(linspace(-1,1,N), linspace(-1,1,N));
signal_amp = exp(-5*(x.^2 + y.^2)); % 模拟幅度
signal_phase = 2*pi*rand(N,N);      % 随机相位 (未知)
% signal = signal_amp .* exp(1i*signal_phase); % 原始复信号

signal = O;

% 傅里叶幅度 (接收端已知)
F_signal = fft2(signal);
F_amp = abs(F_signal);

% 物平面幅度 (发送端已知)
obj_amp = abs(signal);

%% 2. 初始化 GS 算法
% 用物平面幅度和随机相位初始化
f_k = obj_amp .* exp(1i*2*pi*rand(N,N));

max_iter = 20;  % 最大迭代次数
error_history = zeros(max_iter,1);

%% 3. 迭代过程
for k = 1:max_iter
    % 傅里叶变换
    F_k = fft2(f_k);
    
    % 傅里叶幅度约束
    F_k = F_amp .* exp(1i*angle(F_k));
    
    % 逆傅里叶变换
    f_k = ifft2(F_k);
    
    % 物平面幅度约束
    f_k = obj_amp .* exp(1i*angle(f_k));
    
    % 计算误差 (傅里叶幅度)
    error_history(k) = norm(abs(fft2(f_k)) - F_amp, 'fro') / norm(F_amp, 'fro');
end

%% 4. 显示结果
figure;
subplot(2,3,1); imagesc(O_amp); axis image; title('原始幅度'); colorbar;
subplot(2,3,2); imagesc(O_ang); axis image; title('原始相位'); colorbar;
subplot(2,3,3); imagesc(abs(f_k)); axis image; title('恢复幅度'); colorbar;
subplot(2,3,4); imagesc(angle(f_k)); axis image; title('恢复相位'); colorbar;
subplot(2,3,5); plot(1:max_iter,error_history,'LineWidth',2);
xlabel('迭代次数'); ylabel('归一化误差'); title('收敛曲线');

%% 5. 输出相位误差
phase_error = angle(exp(1i*(signal_phase - angle(f_k))));
fprintf('相位均方误差: %.4f radians\n', sqrt(mean(phase_error(:).^2)));