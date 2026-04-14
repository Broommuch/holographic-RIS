function [theta_est, phi_est, alpha_est] = GS_channel_estimation(Z, S, c, N, L, d, maxIter)

% ==========================================================
% GS-based Channel Estimation
% 输入：
%   Z  : N × Tp 强度测量
%   S  : L × Tp 导频
%   c  : N × 1 参考信号
%   N  : RIS单元数
%   L  : 用户数
%   d  : 阵元间距（lambda）
%   maxIter : 最大迭代次数
%
% 输出：
%   theta_est : 方位角
%   phi_est   : 俯仰角（当前未用，可扩展）
%   alpha_est : 复信道增益
% ==========================================================

[N_check, Tp] = size(Z);
assert(N_check == N, 'Dimension mismatch');

%% ================= 初始化 =================
theta_est = rand(L,1)*pi - pi/2;   % 随机初始化
phi_est   = zeros(L,1);            % 先不估计俯仰角（可扩展）
alpha_est = ones(L,1);             % 初始信道

Y = sqrt(Z) .* exp(1j*2*pi*rand(N,Tp));  % 随机相位初始化

%% ================= GS迭代 =================
for iter = 1:maxIter
    
    % ---------- Step 1: 信号投影（LS求解 alpha） ----------
    
    % 构造阵列矩阵 V
    V = zeros(N,L);
    for l = 1:L
        for n = 1:N
            V(n,l) = exp(1j*2*pi*d*(n-1)*sin(theta_est(l)));
        end
    end
    
    % 最小二乘更新 alpha
    for t = 1:Tp
        y_t = Y(:,t) - c;
        alpha_est = (V * S(:,t)) \ y_t;   % LS（简单版）
    end
    
    % ---------- Step 2: 更新 Y ----------
    
    for t = 1:Tp
        y_model = V * (alpha_est .* S(:,t)) + c;
        Y(:,t) = y_model;
    end
    
    % ---------- Step 3: 幅度投影 ----------
    
    Y = sqrt(Z) .* exp(1j * angle(Y));
    
    % ---------- Step 4: 粗角度更新（简单扫描，可替换优化） ----------
    
    theta_grid = linspace(-pi/2, pi/2, 180);
    best_err = inf;
    
    for th = theta_grid
        V_test = zeros(N,L);
        for l = 1:L
            for n = 1:N
                V_test(n,l) = exp(1j*2*pi*d*(n-1)*sin(th));
            end
        end
        
        err = 0;
        for t = 1:Tp
            y_test = V_test * (alpha_est .* S(:,t)) + c;
            err = err + norm(abs(y_test).^2 - Z(:,t))^2;
        end
        
        if err < best_err
            best_err = err;
            theta_est(:) = th;
        end
    end
    
    % 可打印收敛情况
    % fprintf('Iter %d, error = %.4f\n', iter, best_err);
end

end