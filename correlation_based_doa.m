% 算法2: 基于相关性的最大似然估计
function theta_est = correlation_based_doa(P_measured, codebook, angles, lambda, d)
    N_ris = size(codebook, 1);
    correlation = zeros(length(angles), 1);
    
    for i = 1:length(angles)
        theta_test = angles(i);
        a_theta = exp(1j * 2*pi * d/lambda * (0:N_ris-1)' * sind(theta_test));
        
        % 计算理论功率模式
        P_theoretical = zeros(size(P_measured));
        for k = 1:length(P_measured)
            w_k = codebook(:, k);
            P_theoretical(k) = abs(w_k' * a_theta)^2;
        end
        
        % 计算相关性
        correlation(i) = corr(P_measured, P_theoretical);
    end
    
    [~, max_idx] = max(correlation);
    theta_est = angles(max_idx);
end