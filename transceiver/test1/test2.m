% holographic_ris_recovery.m
% 完整仿真：Gauss-Newton / LM 恢复复标量 s (幅度与相位)
clear; close all; clc;

rng(1); % 固定随机种子，便于复现

%% 仿真参数
N = 8; M = 8;         % RIS 尺寸示例 => P = N*M
P = N * M;
SNR_dB = 20;          % 单个观测的信噪比 (dB) for y (signal power vs noise)
numTrials = 20;       % Monte Carlo 重复次数
maxIter = 100;
tol = 1e-8;

% 记录结果
mse_amp = zeros(numTrials,1);
mse_phase = zeros(numTrials,1);
converge_iters = zeros(numTrials,1);

for trial = 1:numTrials
    %% 生成通道 a = phi * h
    phi = 1; % 可改为复常数或未知
    h = (randn(P,1)+1j*randn(P,1))/sqrt(2); % CN(0,1)
    a = phi * h;
    
    %% 参考波 b（已知） -- 示例：随机复数但有空间相位差
    b_amp = 0.5 + rand(P,1)*1.0; % 非零幅度
    b_phase = 2*pi*rand(P,1);
    b = b_amp .* exp(1j*b_phase);
    
    %% 真实发送信号 s_true
    s_amp_true = 1.2; % 真正幅度
    s_phase_true = 0.7; % 真正相位 (rad)
    s_true = s_amp_true * exp(1j * s_phase_true);
    
    %% 生成无噪声量测并加噪声
    y_clean = abs(a * s_true + b).^2;
    % 设定观测噪声功率，使得 SNR_dB 对应 y_clean 的平均功率
    signal_power = mean(y_clean);
    SNR_lin = 10^(SNR_dB/10);
    sigma2 = signal_power / SNR_lin;
    noise = sqrt(sigma2) * randn(P,1); % 观测上为实高斯噪声（幅度平方测量噪声）
    y = y_clean + noise;
    
    %% 初始化 s0（见文中方法A）
    z_hat_mag = sqrt(max(y - sigma2, 0));          % 估计 |z_i|
    z_hat = z_hat_mag .* exp(1j * angle(b));       % 用 b 的相位近似 z 的相位
    denom = sum(abs(a).^2);
    if denom == 0
        s0 = 0;
    else
        s0 = sum(conj(a) .* (z_hat - b)) / denom;  % 线性最小二乘初始化
    end
    
    % 若初始化非常弱，可多次随机化（此处不启用）
    
    %% 调用 Gauss-Newton / LM 迭代求解
    opts.maxIter = maxIter;
    opts.tol = tol;
    opts.lambda0 = 1e-3;
    opts.verbose = false;
    [s_est, info] = gn_lm(y, a, b, s0, opts);
    
    %% 记录结果
    mse_amp(trial) = (abs(s_est) - abs(s_true))^2;
    % 相位误差取主值在 [-pi,pi]
    phase_err = angle(s_est) - angle(s_true);
    phase_err = mod(phase_err+pi, 2*pi) - pi;
    mse_phase(trial) = phase_err^2;
    converge_iters(trial) = info.iter;
    
    if opts.verbose
        fprintf('Trial %d: s_true=%.3f∠%.3f, s_est=%.3f∠%.3f, iter=%d\n', ...
            trial, abs(s_true), angle(s_true), abs(s_est), angle(s_est), info.iter);
    end
end

%% 输出统计
fprintf('--- Results over %d trials ---\n', numTrials);
fprintf('Mean amplitude MSE = %.3e, std = %.3e\n', mean(mse_amp), std(mse_amp));
fprintf('Mean phase MSE (rad^2) = %.3e, std = %.3e\n', mean(mse_phase), std(mse_phase));
fprintf('Mean iterations = %.2f\n', mean(converge_iters));

%% 可视化：最后一次试验的收敛曲线
figure;
plot(info.resnorm_hist, '-o');
xlabel('Iteration'); ylabel('Residual norm (sum r_i^2)');
title('GN/LM Convergence (last trial)'); grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 子函数：Gauss-Newton with Levenberg-Marquardt damping
function [s_est, info] = gn_lm(y, a, b, s0, opts)
% Inputs:
%   y   : P x 1 real measurements
%   a,b : P x 1 complex known vectors
%   s0  : initial complex guess
%   opts: structure with fields:
%         maxIter, tol, lambda0, verbose
% Outputs:
%   s_est: estimated complex scalar
%   info: diagnostics: iter, resnorm_hist

    if ~isfield(opts,'maxIter'), opts.maxIter = 100; end
    if ~isfield(opts,'tol'), opts.tol = 1e-8; end
    if ~isfield(opts,'lambda0'), opts.lambda0 = 1e-3; end
    if ~isfield(opts,'verbose'), opts.verbose = false; end

    P = numel(y);
    s = s0;
    lambda = opts.lambda0;
    resnorm_hist = [];
    prev_resnorm = inf;
    
    for iter = 1:opts.maxIter
        z = a * s + b;                     % P x 1 complex
        z_conj = conj(z);
        % residuals r = y - |z|^2
        r = y - abs(z).^2;                 % P x 1 real
        resnorm = sum(r.^2);
        resnorm_hist(end+1) = resnorm;
        if opts.verbose
            fprintf('Iter %d: resnorm=%.6e, lambda=%.3e\n', iter, resnorm, lambda);
        end
        
        % Compute Jacobian J (P x 2)
        % J(:,1) = dr/ds_r = -2*Re(conj(z).*a)
        % J(:,2) = dr/ds_i =  2*Im(conj(z).*a)
        Az = conj(z) .* a; % P x 1 complex
        J = zeros(P,2);
        J(:,1) = -2 * real(Az);
        J(:,2) =  2 * imag(Az);
        
        % Gauss-Newton step (LM damping)
        H = J' * J;               % 2x2
        g = J' * r;               % 2x1
        % Solve (H + lambda I) delta = g
        A = H + lambda * eye(2);
        delta = A \ g;            % 2x1 real (delta for [s_r; s_i])
        
        % Update candidate
        s_candidate = (real(s) + delta(1)) + 1j*(imag(s) + delta(2));
        
        % Evaluate candidate residual norm
        z_cand = a * s_candidate + b;
        r_cand = y - abs(z_cand).^2;
        resnorm_cand = sum(r_cand.^2);
        
        if resnorm_cand < resnorm
            % Accept step, decrease lambda
            s = s_candidate;
            lambda = lambda / 10;
        else
            % Reject step, increase lambda
            lambda = lambda * 10;
        end
        
        % Check convergence (on parameter step size and residual)
        if norm(delta) < opts.tol * (1 + norm([real(s); imag(s)]))
            break;
        end
        if abs(prev_resnorm - resnorm) < opts.tol * (1 + resnorm)
            break;
        end
        prev_resnorm = resnorm;
    end

    s_est = s;
    info.iter = iter;
    info.resnorm_hist = resnorm_hist;
end
