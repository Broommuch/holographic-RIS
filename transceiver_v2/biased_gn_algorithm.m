function [s_est, info] = biased_gn_algorithm(z_obs, P, b_ref, maxIter)
% BIASED_GN_ALGORITHM
% 用高斯牛顿/LM方法求解带参考项的幅度恢复问题：
%
%       z_obs = abs(P*s + b)
%
% 输入：
%   z_obs   : N_obs x 1，观测幅值，不是能量
%   P       : N_obs x N_sym，已知成形滤波/观测矩阵
%   b_ref   : 已知参考项
%             - 如果 length(b_ref) == N_sym，则认为 b = P*b_ref
%             - 如果 length(b_ref) == N_obs，则认为 b = b_ref
%   maxIter : 最大迭代次数
%
% 输出：
%   s_est : N_sym x 1，恢复的复符号
%   info  : 迭代信息，包括 loss, grad_norm, step_norm

    %% ================= 基本尺寸处理 =================
    z_obs = z_obs(:);
    [N_obs, N_sym] = size(P);

    if length(z_obs) ~= N_obs
        error('z_obs 的长度必须等于 P 的行数。');
    end

    b_ref = b_ref(:);

    if length(b_ref) == N_sym
        b = P * b_ref;
    elseif length(b_ref) == N_obs
        b = b_ref;
    else
        error('b_ref 的长度必须等于 N_sym 或 N_obs。');
    end

    z_obs = max(real(z_obs), 0);

    %% ================= 初始化 =================
    % 用强参考近似做初始化：
    % |P*s+b|^2 ≈ |b|^2 + 2 Re{(P*s) b^*}
    % 这一步不是严格精确，但通常比随机初始化稳定。

    q = z_obs.^2 - abs(b).^2;

    C = diag(conj(b)) * P;

    % q ≈ 2 Re{C*s}
    % s = x + j y
    % Re{C*s} = real(C)*x - imag(C)*y
    M = 2 * [real(C), -imag(C)];

    reg_init = 1e-6 * trace(M' * M) / max(1, size(M, 2));
    xy0 = (M' * M + reg_init * eye(2*N_sym)) \ (M' * q);

    x = xy0(1:N_sym);
    y = xy0(N_sym+1:end);

    s = x + 1j*y;

    % 如果初始化异常，退化为小随机初始化
    if any(~isfinite(s)) || norm(s) < 1e-12
        rng(1);
        s = 0.1 * (randn(N_sym,1) + 1j*randn(N_sym,1));
    end

    %% ================= LM-GN 参数 =================
    lambda = 1e-3;
    lambda_up = 10;
    lambda_down = 0.3;

    minLambda = 1e-12;
    maxLambda = 1e12;

    eps_abs = 1e-12;

    info.loss = zeros(maxIter, 1);
    info.grad_norm = zeros(maxIter, 1);
    info.step_norm = zeros(maxIter, 1);
    info.lambda = zeros(maxIter, 1);

    %% ================= 初始损失 =================
    a = P*s + b;
    r = abs(a) - z_obs;
    loss = 0.5 * norm(r)^2;

    %% ================= 主迭代 =================
    for it = 1:maxIter

        a = P*s + b;
        abs_a = abs(a);
        abs_a_safe = max(abs_a, eps_abs);

        r = abs_a - z_obs;

        % ================= 构造实值 Jacobian =================
        % f_n = |a_n|
        % da = P ds
        %
        % df = Re{ conj(a)/|a| * da }
        %
        % 对 x = real(s):
        % Jx = Re{ conj(a)/|a| * P }
        %
        % 对 y = imag(s):
        % Jy = Re{ conj(a)/|a| * jP }

        phase_factor = conj(a) ./ abs_a_safe;     % N_obs x 1

        JP = diag(phase_factor) * P;

        Jx = real(JP);
        Jy = real(1j * JP);

        J = [Jx, Jy];

        g = J' * r;
        H_gn = J' * J;

        grad_norm = norm(g);

        % ================= LM 步长 =================
        D = diag(diag(H_gn));
        if norm(D, 'fro') < eps_abs
            D = eye(2*N_sym);
        end

        accepted = false;

        for ls = 1:20

            H_lm = H_gn + lambda * D + 1e-12 * eye(2*N_sym);

            delta = - H_lm \ g;

            dx = delta(1:N_sym);
            dy = delta(N_sym+1:end);

            s_trial = s + dx + 1j*dy;

            a_trial = P*s_trial + b;
            r_trial = abs(a_trial) - z_obs;
            loss_trial = 0.5 * norm(r_trial)^2;

            if loss_trial < loss
                s = s_trial;
                loss = loss_trial;
                lambda = max(lambda * lambda_down, minLambda);
                accepted = true;
                break;
            else
                lambda = min(lambda * lambda_up, maxLambda);
            end
        end

        info.loss(it) = loss;
        info.grad_norm(it) = grad_norm;
        info.step_norm(it) = norm(delta);
        info.lambda(it) = lambda;

        % ================= 停止条件 =================
        if ~accepted
            break;
        end

        if norm(delta) / max(norm([real(s); imag(s)]), 1) < 1e-8
            break;
        end

        if grad_norm < 1e-8
            break;
        end
    end

    %% ================= 输出整理 =================
    s_est = s;

    info.loss = info.loss(1:it);
    info.grad_norm = info.grad_norm(1:it);
    info.step_norm = info.step_norm(1:it);
    info.lambda = info.lambda(1:it);
    info.n_iter = it;
    info.final_loss = loss;
    info.final_rel_residual = norm(abs(P*s_est + b) - z_obs) / max(norm(z_obs), eps_abs);

end