%% 1. 系统参数初始化
clear; clc; close all;

% 基本参数
N = 8;              % RIS行数
M = 8;              % RIS列数
NM = N * M;         % RIS总单元数
SNR_dB = 20;        % 信噪比(dB)
num_iter = 100;     % Gauss-Newton迭代次数
num_trials = 50;    % 蒙特卡洛实验次数

% 待恢复信号的真值
s_true = 0.8 + 0.6i;    % 复数信号
s_mag_true = abs(s_true);   % 幅度真值
s_phase_true = angle(s_true); % 相位真值(弧度)

% 反射系数
phi = 1;  % 假设所有单元反射系数相同，可设为1

fprintf('系统参数：\n');
fprintf('  RIS单元数：%d × %d = %d\n', N, M, NM);
fprintf('  信号真值：幅度=%.2f, 相位=%.2f rad (%.1f°)\n', ...
    s_mag_true, s_phase_true, rad2deg(s_phase_true));
fprintf('  SNR：%d dB, 实验次数：%d\n', SNR_dB, num_trials);

%% 2. 信道和参考波生成
% 设置随机种子确保可重复性
rng(2025);

% 生成用户到RIS的信道h (瑞利衰落)
h = (randn(NM, 1) + 1i * randn(NM, 1)) / sqrt(2);

% 生成已知参考波信号b
% 可选方案1：随机QPSK信号
b = exp(1i * 2*pi * rand(NM, 1));

% 可选方案2：线性相位波前
% b = exp(1i * linspace(0, 2*pi, NM)');

% 可选方案3：平面波，波达方向θ_b
% theta_b = pi/6;  % 30度
% array_pos = (0:NM-1)' * 0.5;  % 假设半波长间距
% b = exp(1i * 2*pi * array_pos * sin(theta_b));

fprintf('信道与参考波生成完成。\n');

%% 3. 生成带噪声的观测数据
% 计算无噪声的电场
E_noiseless = phi * h * s_true + b;

% 计算无噪声观测(强度)
y_noiseless = abs(E_noiseless).^2;

% 计算噪声功率
signal_power = mean(y_noiseless);
noise_power = signal_power / (10^(SNR_dB/10));

% 生成高斯白噪声
noise = sqrt(noise_power/2) * (randn(NM, 1) + 1i*randn(NM, 1));
noise = real(noise);  % 只取实部，因为观测是实数

% 添加噪声的观测
y = y_noiseless + noise;

fprintf('观测数据生成完成。\n');
fprintf('  无噪声观测平均功率：%.4f\n', mean(y_noiseless));
fprintf('  噪声功率：%.4f, SNR实际值：%.2f dB\n', ...
    noise_power, 10*log10(signal_power/noise_power));

%% 5. 主仿真：多次实验统计性能
fprintf('\n=== 开始蒙特卡洛仿真 ===\n');

% 存储结果
mag_errors = zeros(num_trials, 1);
phase_errors = zeros(num_trials, 1);
mse_dB = zeros(num_trials, 1);
convergence_flags = zeros(num_trials, 1);
iter_counts = zeros(num_trials, 1);

for trial = 1:num_trials
    fprintf('实验 %d/%d: ', trial, num_trials);
    
    % 重新生成噪声（信道和参考波可固定或变化）
    if mod(trial, 10) == 0
        % 每10次实验重新生成信道，模拟时变
        h = (randn(NM, 1) + 1i * randn(NM, 1)) / sqrt(2);
    end
    
    % 生成带噪声的观测
    noise = sqrt(noise_power/2) * (randn(NM, 1) + 1i*randn(NM, 1));
    y = y_noiseless + real(noise);
    
    % 随机初始化（幅度在0-2，相位在0-2π）
%     mag_init = 2 * rand();
    mag_init = abs(s_true);
%     phase_init = 2 * pi * rand();
    phase_init = angle(s_true);
    s_init = mag_init * exp(1i * phase_init);
    
    % 运行Gauss-Newton算法
    [s_est, s_history, cost_history] = gauss_newton_ris(...
        y, h, b, phi, s_init, num_iter, 1e-6);
    
    % 记录性能
    mag_est = abs(s_est);
    phase_est = angle(s_est);
    
    mag_errors(trial) = abs(mag_est - s_mag_true) / s_mag_true;
    phase_errors(trial) = min(mod(phase_est - s_phase_true, 2*pi), ...
                              mod(s_phase_true - phase_est, 2*pi));
    mse_dB(trial) = 20*log10(abs(s_est - s_true)/abs(s_true));
    
    % 记录收敛信息
    iter_counts(trial) = length(cost_history) - 1;
    convergence_flags(trial) = (iter_counts(trial) < num_iter);
    
    fprintf('估计: 幅度=%.3f(误差%.1f%%), 相位=%.2f rad(误差%.2f rad)\n', ...
        mag_est, 100*mag_errors(trial), phase_est, phase_errors(trial));
end

%% 6. 性能分析和可视化
fprintf('\n=== 性能统计 ===\n');
fprintf('幅度相对误差: 均值=%.2f%%, 标准差=%.2f%%\n', ...
    100*mean(mag_errors), 100*std(mag_errors));
fprintf('相位误差: 均值=%.3f rad, 标准差=%.3f rad\n', ...
    mean(phase_errors), std(phase_errors));
fprintf('MSE(dB): 均值=%.2f dB, 标准差=%.2f dB\n', ...
    mean(mse_dB), std(mse_dB));
fprintf('收敛率: %.1f%% (%d/%d)\n', ...
    100*mean(convergence_flags), sum(convergence_flags), num_trials);
fprintf('平均迭代次数: %.1f\n', mean(iter_counts));

% 绘制收敛曲线（最后一次实验）
figure('Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
plot(0:length(cost_history)-1, 10*log10(cost_history), 'b-', 'LineWidth', 2);
xlabel('迭代次数');
ylabel('代价函数(dB)');
title('代价函数收敛曲线');
grid on;

subplot(1, 3, 2);
s_r_history = s_history(1, :);
s_i_history = s_history(2, :);
plot(s_r_history, s_i_history, 'r.-', 'LineWidth', 1.5, 'MarkerSize', 10);
hold on;
plot(real(s_true), imag(s_true), 'go', 'MarkerSize', 10, 'LineWidth', 2);
plot(real(s_init), imag(s_init), 'ks', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('实部');
ylabel('虚部');
title('估计轨迹');
legend('估计轨迹', '真值', '初始值', 'Location', 'best');
grid on;
axis equal;

subplot(1, 3, 3);
% 绘制幅度误差分布
yyaxis left;
histogram(100*mag_errors, 20, 'FaceColor', 'b', 'FaceAlpha', 0.6);
xlabel('幅度相对误差(%)');
ylabel('频数');
yyaxis right;
histogram(phase_errors, 20, 'FaceColor', 'r', 'FaceAlpha', 0.6);
ylabel('频数');
title('误差分布');
legend('幅度误差', '相位误差');
grid on;

sgtitle('RIS相位检索仿真结果');

% 绘制MSE随SNR变化曲线
figure('Position', [100, 600, 800, 400]);
SNR_range = 0:5:30;
mse_snr = zeros(length(SNR_range), 1);

for i = 1:length(SNR_range)
    snr = SNR_range(i);
    noise_power_snr = signal_power / (10^(snr/10));
    
    % 运行10次取平均
    mse_temp = 0;
    for k = 1:10
        noise = sqrt(noise_power_snr/2) * (randn(NM, 1) + 1i*randn(NM, 1));
        y_snr = y_noiseless + real(noise);
        s_est = gauss_newton_ris(y_snr, h, b, phi, s_init, num_iter, 1e-6);
        mse_temp = mse_temp + 20*log10(abs(s_est - s_true)/abs(s_true));
    end
    mse_snr(i) = mse_temp / 10;
end

plot(SNR_range, mse_snr, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('SNR (dB)');
ylabel('归一化MSE (dB)');
title('MSE vs SNR');
grid on;
hold on;
% 理论下界参考
theoretical = -SNR_range;
plot(SNR_range, theoretical, 'r--', 'LineWidth', 1.5);
legend('仿真MSE', '理论下界(参考)');

%% 4. Gauss-Newton算法实现
function [s_est, s_history, cost_history] = gauss_newton_ris(y, h, b, phi, s_init, max_iter, tol)
    % Gauss-Newton法求解RIS相位检索问题
    %
    % 输入：
    %   y: 观测向量 (NM×1)
    %   h: 信道向量 (NM×1)
    %   b: 参考波向量 (NM×1)
    %   phi: 反射系数标量
    %   s_init: 初始估计 (复数)
    %   max_iter: 最大迭代次数
    %   tol: 收敛容差
    %
    % 输出：
    %   s_est: 估计的信号
    %   s_history: 迭代历史
    %   cost_history: 代价函数历史
    
    % 将复数s分解为实部虚部
    x = [real(s_init); imag(s_init)];  % 待优化变量 [s_r; s_i]
    
    % 预计算常数
    NM = length(y);
    A = phi * h;  % 预处理
    
    s_history = zeros(2, max_iter+1);
    s_history(:, 1) = x;
    cost_history = zeros(max_iter+1, 1);
    
    % 计算初始代价
    residuals = compute_residuals(x, y, A, b);
    cost_history(1) = sum(residuals.^2);
    
    fprintf('Gauss-Newton迭代开始，初始代价：%.6e\n', cost_history(1));
    
    for iter = 1:max_iter
        % 计算当前残差
        residuals = compute_residuals(x, y, A, b);
        
        % 计算雅可比矩阵
        J = compute_jacobian(x, A, b);
        
        % Gauss-Newton更新：Δx = -(J'J)^{-1} J' r
        grad = J' * residuals;  % 梯度
        Hessian_approx = J' * J;  % 近似Hessian
        
        % 添加正则化防止奇异
        mu = 1e-6 * eye(2);
        delta_x = - (Hessian_approx + mu) \ grad;
        
        % 更新估计
        x_new = x + delta_x;
        
        % 计算新代价
        residuals_new = compute_residuals(x_new, y, A, b);
        cost_new = sum(residuals_new.^2);
        
        % 检查是否接受更新（简单线搜索）
        alpha = 1.0;
        while cost_new > cost_history(iter) && alpha > 1e-4
            alpha = alpha * 0.5;
            x_new = x + alpha * delta_x;
            residuals_new = compute_residuals(x_new, y, A, b);
            cost_new = sum(residuals_new.^2);
        end
        
        % 更新变量
        x = x_new;
        s_history(:, iter+1) = x;
        cost_history(iter+1) = cost_new;
        
        % 检查收敛
        rel_change = norm(delta_x) / (norm(x) + eps);
        if rel_change < tol
            fprintf('  迭代 %d: 代价=%.6e, 相对变化=%.3e (已收敛)\n', ...
                iter, cost_new, rel_change);
            break;
        end
        
        if mod(iter, 10) == 0
            fprintf('  迭代 %d: 代价=%.6e, 相对变化=%.3e\n', ...
                iter, cost_new, rel_change);
        end
    end
    
    % 转换为复数输出
    s_est = x(1) + 1i * x(2);
    
    % 截断历史记录
    s_history = s_history(:, 1:iter+1);
    cost_history = cost_history(1:iter+1);
end

%% 辅助函数：计算残差
function residuals = compute_residuals(x, y, A, b)
    % 计算残差 r_i = y_i - |A_i * s + b_i|^2
    s = x(1) + 1i * x(2);
    E = A * s + b;
    y_pred = abs(E).^2;
    residuals = y - y_pred;
end

%% 辅助函数：计算雅可比矩阵
function J = compute_jacobian(x, A, b)
    % 计算雅可比矩阵，每行对应一个观测
    s = x(1) + 1i * x(2);
    E = A * s + b;
    
    % 计算u_i和v_i
    u = real(E);
    v = imag(E);
    
    % 计算偏导数
    du_dsr = real(A);
    dv_dsr = imag(A);
    du_dsi = -imag(A);  % real(j*A) = -imag(A)
    dv_dsi = real(A);   % imag(j*A) = real(A)
    
    % 计算雅可比矩阵行：J_i = [dr_i/ds_r, dr_i/ds_i]
    NM = length(A);
    J = zeros(NM, 2);
    
    for i = 1:NM
        J(i, 1) = -2 * (u(i) * du_dsr(i) + v(i) * dv_dsr(i));
        J(i, 2) = -2 * (u(i) * du_dsi(i) + v(i) * dv_dsi(i));
    end
end

%% 7. 扩展：Levenberg-Marquardt算法（更鲁棒）
function [s_est, lambda_history, cost_history] = levenberg_marquardt_ris(y, h, b, phi, s_init, max_iter, tol)
    % Levenberg-Marquardt算法
    x = [real(s_init); imag(s_init)];
    NM = length(y);
    A = phi * h;
    
    lambda = 0.01;  % 阻尼因子初始值
    lambda_history = zeros(max_iter, 1);
    cost_history = zeros(max_iter+1, 1);
    
    % 初始代价
    residuals = compute_residuals(x, y, A, b);
    cost_history(1) = sum(residuals.^2);
    
    for iter = 1:max_iter
        % 计算雅可比和梯度
        J = compute_jacobian(x, A, b);
        grad = J' * residuals;
        H_approx = J' * J;
        
        % Levenberg-Marquardt更新
        delta_x = - (H_approx + lambda * diag(diag(H_approx))) \ grad;
        
        % 试探更新
        x_new = x + delta_x;
        residuals_new = compute_residuals(x_new, y, A, b);
        cost_new = sum(residuals_new.^2);
        
        % 判断是否接受更新
        rho = (cost_history(iter) - cost_new) / ...
              (delta_x' * (lambda * diag(diag(H_approx)) * delta_x + grad));
        
        if rho > 0
            % 接受更新
            x = x_new;
            cost_history(iter+1) = cost_new;
            lambda = lambda * max(1/3, 1 - (2*rho - 1)^3);
        else
            % 拒绝更新
            cost_history(iter+1) = cost_history(iter);
            lambda = lambda * 2;
        end
        
        lambda_history(iter) = lambda;
        
        % 检查收敛
        if norm(delta_x) < tol * (norm(x) + eps)
            break;
        end
    end
    
    s_est = x(1) + 1i * x(2);
    cost_history = cost_history(1:iter+1);
    lambda_history = lambda_history(1:iter);
end