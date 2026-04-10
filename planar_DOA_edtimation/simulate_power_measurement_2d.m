%% 4. 功率测量模拟
function P_measured = simulate_power_measurement_2d(codebook, azimuth, elevation, ...
                                                  M, N, d, lambda, SNR_dB)
% 模拟二维阵列的功率测量
    
    K = size(codebook, 2);
    P_measured = zeros(K, 1);
    
    % 生成导向矢量
    a_true = steering_vector_2d(azimuth, elevation, M, N, d, lambda);
    
    for k = 1:K
        w_k = codebook(:, k);
        
        % 接收信号
        y_k = w_k' * a_true;
        
        % 添加噪声
        noise_power = 10^(-SNR_dB/10) * abs(y_k)^2;
        y_k_noisy = y_k + sqrt(noise_power/2) * (randn(1) + 1j*randn(1));
        
        P_measured(k) = abs(y_k_noisy)^2;
    end
end