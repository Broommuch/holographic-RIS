function [s_best, loss_best, loss_hist] = biased_gs_qpsk(z, A, b, const, n_iter, n_restart)
% 模型: z = abs(A*s + b)
% A: N_obs x N_sym
% b: N_obs x 1
% z: N_obs x 1
% const: QPSK星座列向量

    [N_obs, N_sym] = size(A);

    loss_best = inf;
    s_best = zeros(N_sym, 1);
    loss_hist = [];

    for rr = 1:n_restart

        % 随机QPSK初始化
        idx0 = randi(length(const), N_sym, 1);
        s = const(idx0);

        for it = 1:n_iter

            u = A * s + b;

            % 避免除0
            u_abs = abs(u);
            u_abs(u_abs < 1e-12) = 1e-12;

            % 幅度投影
            u_proj = z .* u ./ u_abs;

            % 反投影到符号域：min ||A*s + b - u_proj||^2
            s_ls = A \ (u_proj - b);

            % QPSK硬投影
            s_new = zeros(N_sym, 1);
            for k = 1:N_sym
                [~, id] = min(abs(s_ls(k) - const));
                s_new(k) = const(id);
            end

            s = s_new;

            loss = norm(abs(A*s + b) - z)^2 / norm(z)^2;

            if loss < loss_best
                loss_best = loss;
                s_best = s;
            end
        end
    end
end