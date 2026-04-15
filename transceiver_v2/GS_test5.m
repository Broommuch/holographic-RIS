%% 数字全息仿真示例（带 GS 迭代）依旧错误，没用
clc; clear; close all;

%% 1. 生成目标信号 (平面波入射)
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

%% 2. 生成参考平面波
lambda = 633e-9;
pixel_size = 5e-6;

theta_x_ref = 10; % x方向倾斜角
theta_y_ref = 5;  % y方向倾斜角
theta_x_rad = deg2rad(theta_x_ref);
theta_y_rad = deg2rad(theta_y_ref);

[X_idx, Y_idx] = meshgrid(0:N-1, 0:N-1);
X = X_idx * pixel_size;
Y = Y_idx * pixel_size;

R = exp(1i*2*pi/lambda * (X*sin(theta_x_rad) + Y*sin(theta_y_rad)));

%% 3. 模拟干涉强度
I = abs(O + R).^2;

%% 4. 显示原始目标和干涉图
figure;
subplot(2,3,1); imagesc(O_amp); axis image; colorbar; title('目标幅度 |O|');
subplot(2,3,2); imagesc(angle(O)); axis image; colorbar; title('目标相位 φ_O');
subplot(2,3,3); imagesc(I); axis image; colorbar; title('干涉强度 I');

%% 5. 初始化 GS 迭代
num_iter = 50; % 迭代次数
% O_est = ones(N,N) .* exp(1i*2*pi*rand(N,N)); % 随机初相位

O_est = angle(O);

%% 6. GS 迭代
for iter = 1:num_iter
    % 前向：生成全息
    H = O_est + R;
    
    % 幅度约束：使用测量强度
    H_new = sqrt(I) .* exp(1i*angle(H));
    
    % 回到目标平面
    O_est = H_new - R;
    
    % 可选：每10步显示恢复情况
    if mod(iter,10)==0 || iter==1
        fprintf('迭代 %d / %d\n', iter, num_iter);
        subplot(2,3,5); imagesc(abs(O_est)); axis image; colorbar; title(['恢复幅度 |O|_{rec} (iter=' num2str(iter) ')']);
        subplot(2,3,6); imagesc(angle(O_est)); axis image; colorbar; title(['恢复相位 φ_{rec} (iter=' num2str(iter) ')']);
        drawnow;
    end
end

%% 7. 最终恢复显示
subplot(2,3,5); imagesc(abs(O_est)); axis image; colorbar; title(['最终恢复幅度 |O|_{rec}']);
subplot(2,3,6); imagesc(angle(O_est)); axis image; colorbar; title(['最终恢复相位 φ_{rec}']);

%%
ang_result = angle(O_est);
true_ang = angle(O);
residual = ang_result - true_ang;
