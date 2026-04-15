%% 这个算法看来也是可行的，但是还是对初值很敏感，也就是说还是在解一个欠定问题，看能不能提取一下结构
clc; clear; close all;

%% ================= 1. 参数 =================
N = 128;
d = 0.5;

theta_true = 30; 
phi_true   = 20;

theta_rad = deg2rad(theta_true);
phi_rad   = deg2rad(phi_true);

u_true = sin(theta_rad)*cos(phi_rad);
v_true = sin(theta_rad)*sin(phi_rad);

[x_idx, y_idx] = meshgrid(0:N-1, 0:N-1);
cx = (N-1)/2; cy = (N-1)/2;

%% ================= 2. 构造目标场 =================
phase_true = 2*pi*d*((x_idx-cx)*u_true + (y_idx-cy)*v_true);
O_true = exp(1i * phase_true);

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

%% ================= 5. GN初始化 =================
u = 0.465;   % 初值（可以离真值不太远）
v = 0.178;

num_iter = 20;

figure;

%% ================= 6. Gauss-Newton =================
for iter = 1:num_iter
    
    % 当前O
    phase = 2*pi*d*((x_idx-cx)*u + (y_idx-cy)*v);
    O = exp(1i * phase);
    
    Z = O + R;
    Z_conj = conj(Z);
    
    % 残差
    r = abs(Z).^2 - I;
    
    % ===== Jacobian =====
    dO_du = 1i * 2*pi*d*(x_idx-cx) .* O;
    dO_dv = 1i * 2*pi*d*(y_idx-cy) .* O;
    
    % d|Z|^2/du
    J_u = 2 * real(dO_du .* Z_conj);
    J_v = 2 * real(dO_dv .* Z_conj);
    
    % reshape成向量
    J = [J_u(:), J_v(:)];   % (N^2 × 2)
    r_vec = r(:);
    
    % ===== Gauss-Newton更新 =====
    delta = -(J' * J) \ (J' * r_vec);
    
    u = u + delta(1);
    v = v + delta(2);
    
    % ===== 显示 =====
    fprintf('Iter %d: u=%.4f, v=%.4f\n', iter, u, v);
    
    subplot(1,2,1);
    imagesc(angle(O)); axis image; colorbar;
    title(['恢复相位 iter=' num2str(iter)]);
    
    subplot(1,2,2);
    phase_err = angle(O .* conj(O_true));
    imagesc(phase_err); axis image; colorbar;
    title('相位误差');
    
    drawnow;
end

%% ================= 7. 结果对比 =================
fprintf('\n真实: u=%.4f, v=%.4f\n', u_true, v_true);
fprintf('估计: u=%.4f, v=%.4f\n', u, v);

figure;
imagesc(angle(O .* conj(O_true)));
axis image; colorbar;
title('最终相位误差');