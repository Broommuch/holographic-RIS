% 算法3: 基于压缩感知的稀疏重构
function theta_est = compressive_sensing_doa(P_measured, codebook, angles, lambda, d)
    N_ris = size(codebook, 1);
    L = length(angles);
    
    % 构建感知矩阵
    A = zeros(length(P_measured), L);
    for i = 1:L
        a_theta = exp(1j * 2*pi * d/lambda * (0:N_ris-1)' * sind(angles(i)));
        for k = 1:length(P_measured)
            w_k = codebook(:, k);
            A(k, i) = abs(w_k' * a_theta)^2;
        end
    end
    
    % L1最小化稀疏重构
    cvx_begin quiet
        variable x(L)
        minimize( norm(x, 1) )
        subject to
            norm(A * x - P_measured, 2) <= 0.1 * norm(P_measured, 2)
            x >= 0
    cvx_end
    
    [~, max_idx] = max(x);
    theta_est = angles(max_idx);
end