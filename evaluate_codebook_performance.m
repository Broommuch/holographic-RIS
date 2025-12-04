%% 3. 性能评估函数
function [theta_est, error, time_elapsed] = evaluate_codebook_performance(...
    codebook, theta_true, angles, SNR_dB, lambda, d)
% 评估特定码本的性能
    
    N_ris = size(codebook, 1);
    K = size(codebook, 2);
    
    % 导向矢量
    steering_vector = @(theta) exp(1j * 2*pi * d/lambda * (0:N_ris-1)' * sind(theta));
    s_true = steering_vector(theta_true);
    
    % 功率测量（带噪声）
    P_measured = zeros(K, 1);
    for k = 1:K
        w_k = codebook(:, k);
        y_k = w_k' * s_true;
        noise_power = 10^(-SNR_dB/10) * abs(y_k)^2;
        y_k_noisy = y_k + sqrt(noise_power/2) * (randn(1) + 1j*randn(1));
        P_measured(k) = abs(y_k_noisy)^2;
    end
    
    % 相关性DOA估计
    tic;
    correlation = zeros(length(angles), 1);
    for i = 1:length(angles)
        theta_test = angles(i);
        a_theta = steering_vector(theta_test);
        
        P_theoretical = zeros(K, 1);
        for k = 1:K
            w_k = codebook(:, k);
            P_theoretical(k) = abs(w_k' * a_theta)^2;
        end
        
        correlation(i) = corr(P_measured, P_theoretical);
    end
    
    [~, max_idx] = max(correlation);
    theta_est = angles(max_idx);
    time_elapsed = toc;
    error = abs(theta_est - theta_true);
end