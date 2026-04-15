%% 这个算法初步验证了wirtinger flow算法的可行性，但是不适用当前问题
clc; clear; close all;

%% ================= 1. 参数 =================
N = 128;
d = 0.5;

theta = 30; phi = 20;
theta_rad = deg2rad(theta);
phi_rad   = deg2rad(phi);

u = sin(theta_rad)*cos(phi_rad);
v = sin(theta_rad)*sin(phi_rad);

[x_idx, y_idx] = meshgrid(0:N-1, 0:N-1);
cx = (N-1)/2; cy = (N-1)/2;

%% ================= 2. 目标场 =================
phase_O = 2*pi*d*((x_idx-cx)*u + (y_idx-cy)*v);
O_true = exp(1i * phase_O);

%% ================= 3. 参考光 =================
lambda = 633e-9;
pixel_size = 5e-6;

theta_x_ref = 20;
theta_y_ref = 10;

[X_idx, Y_idx] = meshgrid(0:N-1, 0:N-1);
X = (X_idx-cx)*pixel_size;
Y = (Y_idx-cy)*pixel_size;

R = exp(1i*2*pi/lambda*(X*sin(deg2rad(theta_x_ref)) + ...
                        Y*sin(deg2rad(theta_y_ref))));

%% ================= 4. 测量 =================
I = abs(O_true + R).^2;

%% ================= 5. WF初始化 =================
rng(1);
X = exp(1i * 2*pi*rand(N,N));  % 随机初始化

num_iter = 500;
mu = 0.1;   % 步长（可调）

loss_curve = zeros(num_iter,1);

figure;

%% ================= 6. Wirtinger Flow =================
for iter = 1:num_iter
    
    % 当前估计
    Z = X + R;
    
    % 误差
    residual = (abs(Z).^2 - I);
    
    % Wirtinger梯度
    grad = residual .* Z;
    
    % 梯度下降
    X = X - mu * grad;
    
    % 可选：幅度归一（增强稳定性）
    X = exp(1i * angle(X));
    
    % 记录loss
    loss = 0.5 * sum(residual(:).^2);
    loss_curve(iter) = loss;
    
    % 可视化
    if mod(iter,20)==0 || iter==1
        fprintf('Iter %d, Loss = %.4e\n', iter, loss);
        
        subplot(2,2,1);
        imagesc(angle(X)); axis image; colorbar;
        title(['恢复相位 iter=' num2str(iter)]);
        
        subplot(2,2,2);
        phase_err = angle(X .* conj(O_true));
        imagesc(phase_err); axis image; colorbar;
        title('相位误差');
        
        subplot(2,2,3);
        plot(loss_curve(1:iter),'LineWidth',1.5);
        title('Loss收敛曲线'); grid on;
        
        drawnow;
    end
end

%% ================= 7. 最终结果 =================
figure;
subplot(1,2,1);
imagesc(angle(O_true)); axis image; colorbar;
title('真实相位');

subplot(1,2,2);
imagesc(angle(X)); axis image; colorbar;
title('WF恢复相位');

figure;
imagesc(angle(X .* conj(O_true))); axis image; colorbar;
title('最终相位误差');