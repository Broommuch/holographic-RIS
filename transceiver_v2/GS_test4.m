%% 数字全息仿真示例（改进版）这个版本稍微有点效果，可以作为参考 但是缺少迭代的步骤
clc; clear; close all;

%% 1. 生成目标信号 (平面波入射)
rng(1);

% 参数设置
N = 128;           % 阵列大小
d = 0.5;           % 阵元间距（以波长为单位）
theta = 30;        % 方位角 (度)
phi = 20;          % 俯仰角 (度)

% 转换为弧度
theta_rad = deg2rad(theta);
phi_rad = deg2rad(phi);

% 阵列索引
[x_idx, y_idx] = meshgrid(0:N-1, 0:N-1);

% 中心点
center_x = (N-1)/2;
center_y = (N-1)/2;

% 方向余弦
u = sin(theta_rad) * cos(phi_rad);
v = sin(theta_rad) * sin(phi_rad);

% 相位矩阵
phase_relative = 2*pi * d * ((x_idx - center_x) * u + (y_idx - center_y) * v);

% 目标信号幅度
O_amp = ones(N,N) + (2*rand(N)-1)/10;  % 加一点随机噪声

% 复信号
O = O_amp .* exp(1i * phase_relative);

%% 2. 生成参考平面波
lambda = 633e-9;       % 波长
pixel_size = 5e-6;     % 像素尺寸

% 参考波倾斜角
theta_x_ref = 10; % x方向倾斜角 [度]
theta_y_ref = 5;  % y方向倾斜角 [度]

% 转弧度
theta_x_rad = deg2rad(theta_x_ref);
theta_y_rad = deg2rad(theta_y_ref);

% 像素坐标（物理单位）
[X_idx, Y_idx] = meshgrid(0:N-1, 0:N-1);
X = X_idx * pixel_size;
Y = Y_idx * pixel_size;

% 参考波
R = exp(1i * 2*pi / lambda * (X*sin(theta_x_rad) + Y*sin(theta_y_rad)));

% 可视化参考波
figure;
subplot(1,2,1);
imagesc(real(R)); axis image; colorbar; title('参考波实部');
subplot(1,2,2);
imagesc(angle(R)); axis image; colorbar; colormap hsv; title('参考波相位');

%% 3. 模拟干涉强度
I = abs(O + R).^2;

%% 4. 显示干涉条纹
figure;
subplot(2,3,1); imagesc(O_amp); axis image; colorbar; title('目标幅度 |O|');
subplot(2,3,2); imagesc(angle(O)); axis image; colorbar; title('目标相位 φ_O');
subplot(2,3,3); imagesc(I); axis image; colorbar; title('干涉强度 I');

%% 5. 傅里叶变换 & 频域滤波
F_I = fftshift(fft2(I));

% 归一化频率
kx_norm = sin(theta_x_rad) * lambda / pixel_size; 
ky_norm = sin(theta_y_rad) * lambda / pixel_size;

% FFT索引中心
center = floor(N/2)+1;
cx_plus = round(center + kx_norm * N);
cy_plus = round(center + ky_norm * N);

% 高斯滤波窗口
sigma = 8;
[X, Y] = meshgrid(1:N, 1:N);
Gauss_filter = exp(-((X - cx_plus).^2 + (Y - cy_plus).^2)/(2*sigma^2));

% 频域滤波
F_filtered = F_I .* Gauss_filter;

% 可视化滤波后的频谱
subplot(2,3,4);
imagesc(log(abs(F_filtered)+1)); axis image; colorbar;
title(['高斯滤波后频谱 (中心: ' num2str(cx_plus) ',' num2str(cy_plus) ')']);

% 移频到中心
offset_x = center - cx_plus;
offset_y = center - cy_plus;
F_centered = circshift(F_filtered, [offset_y, offset_x]);

%% 6. 逆傅里叶变换提取复信号
O_recon = ifft2(ifftshift(F_centered));

% 不必再除参考波，已通过频域滤波提取+1级
%% 7. 显示恢复结果
subplot(2,3,5); imagesc(abs(O_recon)); axis image; colorbar; title('恢复幅度 |O|_{rec}');
subplot(2,3,6); imagesc(angle(O_recon)); axis image; colorbar; title('恢复相位 φ_{rec}');