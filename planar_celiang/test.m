%% 压缩波束成形：UPA 单信源 DOA 估计（允许 off-grid）
clear; clc; close all;

% -----------------------------
% 1. 阵列参数
% -----------------------------
Nx = 8;          % x 方向阵元数
Ny = 8;          % y 方向阵元数
N = Nx * Ny;     % 总阵元数

% -----------------------------
% 2. 真实信源参数（off-grid 示例）
% -----------------------------
theta_true = 20 * pi/180;   % 俯仰角 (rad)
phi_true   = -35 * pi/180;  % 方位角 (rad)

% 构建真实导向矢量 a_true = ay ⊗ ax
ax = exp(-1j * pi * (0:Nx-1).' * sin(theta_true) * cos(phi_true));
ay = exp(-1j * pi * (0:Ny-1).' * sin(theta_true) * sin(phi_true));
a_true = kron(ay, ax);  % 注意：kron(ay, ax) 对应 vec(A) = ay ⊗ ax

% -----------------------------
% 3. 构建过完备字典 D (高分辨率角度网格)
% -----------------------------
L_theta = 60;   % 俯仰角网格点数
L_phi   = 60;   % 方位角网格点数
L = L_theta * L_phi;

theta_grid = linspace(0, pi/2, L_theta);      % 0 ~ 90°
phi_grid   = linspace(-pi, pi, L_phi);        % -180° ~ 180°

D = zeros(N, L);  % 字典矩阵
col = 1;
for i = 1:L_theta
    for j = 1:L_phi
        th = theta_grid(i);
        ph = phi_grid(j);
        ax_d = exp(-1j * pi * (0:Nx-1).' * sin(th) * cos(ph));
        ay_d = exp(-1j * pi * (0:Ny-1).' * sin(th) * sin(ph));
        D(:, col) = kron(ay_d, ax_d);
        col = col + 1;
    end
end

% -----------------------------
% 4. 设计测量矩阵 B：随机选取 K 行 2D-DFT 矩阵
% -----------------------------
K = 12;  % 测量次数 (建议 K >= 4, 越大越鲁棒)
snr_db = 20;  % 信噪比 (dB)

% 生成 2D-DFT 矩阵 (注意：MATLAB fft 是按列，需调整)
Fx = dftmtx(Nx) / sqrt(Nx);
Fy = dftmtx(Ny) / sqrt(Ny);
F2D = kron(Fy, Fx);  % 2D-DFT 矩阵 (N x N)

% 随机选择 K 行作为测量向量 (共轭转置用于内积)
rows = randperm(N, K);
B = F2D(rows, :);    % 每行是一个 b^{(k)H}

% -----------------------------
% 5. 模拟测量 y = B * a_true + noise
% -----------------------------
y_clean = B * a_true;
noise_power = norm(y_clean)^2 / (K * 10^(snr_db/10));
w = sqrt(noise_power/2) * (randn(K,1) + 1j*randn(K,1));
y = y_clean + w;

% -----------------------------
% 6. 压缩感知恢复：OMP (1-sparse)
% -----------------------------
A = B * D;  % 感知矩阵 (K x L)

% 归一化列（OMP 要求）
for i = 1:L
    A(:,i) = A(:,i) / norm(A(:,i));
end

% OMP for 1-sparse
residual = y;
max_iter = 1;
support = [];

for iter = 1:max_iter
    % 找与残差最相关的原子
    corr = abs(A' * residual);
    [~, idx] = max(corr);
    support = union(support, idx);
    
    % 最小二乘更新
    x_ls = A(:, support) \ y;
    residual = y - A(:, support) * x_ls;
end

s_est = zeros(L,1);
s_est(support) = x_ls;

% -----------------------------
% 7. 从稀疏解提取角度估计
% -----------------------------
[~, peak_idx] = max(abs(s_est));

% 将线性索引转为 (i,j)
[i_est, j_est] = ind2sub([L_theta, L_phi], peak_idx);

theta_est = theta_grid(i_est);
phi_est   = phi_grid(j_est);

% -----------------------------
% 8. 结果显示
% -----------------------------
fprintf('真实角度: theta = %.2f°, phi = %.2f°\n', ...
        theta_true*180/pi, phi_true*180/pi);
fprintf('估计角度: theta = %.2f°, phi = %.2f°\n', ...
        theta_est*180/pi, phi_est*180/pi);
fprintf('角度误差: d_theta = %.2f°, d_phi = %.2f°\n', ...
        abs(theta_est-theta_true)*180/pi, ...
        min(abs(phi_est-phi_true), 2*pi-abs(phi_est-phi_true))*180/pi);

% 可视化稀疏谱（可选）
figure;
S = reshape(abs(s_est), L_theta, L_phi);
imagesc(phi_grid*180/pi, theta_grid*180/pi, S);
xlabel('Azimuth \phi (deg)'); ylabel('Elevation \theta (deg)');
title('Spatial Spectrum (Compressive Beamforming)');
colorbar;
hold on;
plot(phi_true*180/pi, theta_true*180/pi, 'rx', 'MarkerSize',12, 'LineWidth',2);
plot(phi_est*180/pi,   theta_est*180/pi,   'go', 'MarkerSize',8,  'LineWidth',2);
legend('True', 'Estimated');