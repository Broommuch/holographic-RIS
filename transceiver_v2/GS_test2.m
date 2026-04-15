%% 数字全息仿真示例
clc; clear; close all;

%% 1. 生成目标信号 (模拟通信信号)，平面波入射模型
rng(1);
% 假设N×N均匀平面阵列，阵元间距d（通常为半波长）
% 入射波方向由方位角θ和俯仰角φ定义

% 参数设置
N = 128;           % 阵列大小
d = 0.5;          % 阵元间距（以波长为单位，d=λ/2）
theta = 30;       % 方位角（度），范围[-90, 90]
phi = 20;         % 俯仰角（度），范围[0, 180]

% 转换为弧度
theta_rad = deg2rad(theta);
phi_rad = deg2rad(phi);

% 生成阵列索引
[x_idx, y_idx] = meshgrid(0:N-1, 0:N-1);

% 计算方向余弦
u = sin(theta_rad) * cos(phi_rad);  % x方向余弦
v = sin(theta_rad) * sin(phi_rad);  % y方向余弦

% 计算相对相位延迟
% 相对于阵列中心(N/2, N/2)的相位
center_x = (N-1)/2;
center_y = (N-1)/2;
phase_relative = 2*pi * d * ((x_idx - center_x) * u + (y_idx - center_y) * v);

% 生成目标信号相位矩阵
O_phase = mod(phase_relative, 2*pi);  % 取模到[0, 2π)
% 或直接使用相位值，不需取模
% O_phase = phase_relative;

% 目标信号幅度 (高斯噪声)
O_amp = ones(N,N) + (2*rand(N)-1)/10;

% 复信号
O = O_amp .* exp(1i*O_phase);

%% 2. 生成参考平面波
lambda = 633e-9;
d = 5e-6;  % 像素尺寸

% 归一化频率（在傅里叶变换中范围[-0.5, 0.5]）
% 这对应平面波的倾斜角
kx_norm = 0.2;  % 归一化频率，无量纲
kx = kx_norm;
ky_norm = 0.2;  % 归一化频率，无量纲
ky = kx_norm;

% 从归一化频率计算实际倾斜角
theta_x_actual = asind(kx_norm * lambda / d);  % 实际角度 [度]
theta_y_actual = asind(ky_norm * lambda / d);
fprintf('实际倾斜角: θ_x=%.2f°, θ_y=%.2f°\n', theta_x_actual, theta_y_actual);

% 像素坐标（无量纲索引）
[X_idx, Y_idx] = meshgrid(0:N-1, 0:N-1);

% 参考波（使用归一化频率）
R = exp(1i * 2*pi * (kx_norm * X_idx + ky_norm * Y_idx));

% 可视化
figure;
subplot(1,2,1);
imagesc(real(R));
title('参考波实部');
axis image; colorbar;

subplot(1,2,2);
imagesc(angle(R));
title('参考波相位');
axis image; colorbar; colormap hsv;

%% 3. 模拟干涉强度测量
I = abs(O + R).^2;

%% 4. 显示干涉条纹
figure;
subplot(2,3,1); imagesc(O_amp); axis image; colorbar; title('目标幅度 |O|');
subplot(2,3,2); imagesc(O_phase); axis image; colorbar; title('目标相位 φ_O');
subplot(2,3,3); imagesc(I); axis image; colorbar; title('干涉强度 I');

%% 5. 傅里叶变换 & 频域滤波 (优化版)
F_I = fftshift(fft2(I));

% --- 1. 定义滤波窗口参数 ---
center = floor(N/2) + 1; % 确保中心坐标是整数
% 假设 kx, ky 是归一化的频率偏移 (范围 -0.5 到 0.5)
% 计算 +1 级频谱的中心坐标
cx_plus = center + round(kx * N);
cy_plus = center + round(ky * N);

% 方法1：生成与频域矩阵相同大小的坐标网格
sigma = 8; % 根据你的频谱扩散程度调整，越大包含信息越多，但噪声也越多
center_y = floor(N/2) + 1;
center_x = floor(N/2) + 1;
[X, Y] = meshgrid(1:N, 1:N);  % 生成N×N的网格
Gauss_filter = exp(-((X - cx_plus).^2 + (Y - cy_plus).^2) / (2*sigma^2));

% 应用滤波
F_filtered = F_I .* Gauss_filter;

% --- 3. 显示滤波后的频谱 (验证是否截取正确) ---
subplot(2,3,4); 
imagesc(log(abs(F_filtered)+1)); 
axis image; colorbar; 
title(['高斯滤波后频谱 (中心: ' num2str(cx_plus) ',' num2str(cy_plus) ')']);


% 5. 可选：移频到中心（消除载频）
offset_y = center_y - cy_plus;
offset_x = center_x - cx_plus;
F_centered = circshift(F_filtered, [offset_y, offset_x]);

%% 6. 逆傅里叶变换提取复信号
O_recon = ifft2(ifftshift(F_centered));

% 除以已知参考波
O_recon = O_recon ./ conj(R);

%% 7. 显示恢复结果
subplot(2,3,5); imagesc(abs(O_recon)); axis image; colorbar; title('恢复幅度 |O|_{rec}');
subplot(2,3,6); imagesc(angle(O_recon)); axis image; colorbar; title('恢复相位 φ_{rec}');