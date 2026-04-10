%% 8. 性能分析函数
function analyze_performance_2d(codebook_2d,M,N,d,lambda,azimuth_range,...
    elevation_range,total_elements,K)
% 二维性能分析
    
    fprintf('\n=== 二维性能分析 ===\n');
    
    % 测试不同信噪比
    SNR_test_range = [0, 10, 20, 30];
    n_trials = 50;
    
    results = zeros(length(SNR_test_range), 4); % 存储平均误差
    
    for snr_idx = 1:length(SNR_test_range)
        SNR_test = SNR_test_range(snr_idx);
        errors_az = zeros(n_trials, 1);
        errors_el = zeros(n_trials, 1);
        
        for trial = 1:n_trials
            % 随机生成测试角度
            az_test = -60 + 120*rand();
            el_test = 60*rand();
            
            % 功率测量
            P_test = simulate_power_measurement_2d(codebook_2d, az_test, ...
                                                  el_test, M, N, d, lambda, SNR_test);
            
            % DOA估计
            [az_est, el_est, ~] = correlation_doa_2d(P_test, codebook_2d, ...
                                                    azimuth_range, elevation_range, M, N, d, lambda);
            
            errors_az(trial) = abs(az_est - az_test);
            errors_el(trial) = abs(el_est - el_test);
        end
        
        results(snr_idx, :) = [SNR_test, mean(errors_az), mean(errors_el), ...
                              mean(sqrt(errors_az.^2 + errors_el.^2))];
    end
    
    % 显示性能结果
    figure('Position', [200, 200, 800, 600]);
    subplot(2,2,1);
    semilogy(results(:,1), results(:,2), 'ro-', 'LineWidth', 2);
    hold on;
    semilogy(results(:,1), results(:,3), 'bs-', 'LineWidth', 2);
    semilogy(results(:,1), results(:,4), 'g^-', 'LineWidth', 2);
    xlabel('信噪比 (dB)');
    ylabel('估计误差 (度)');
    title('不同信噪比下的估计误差');
    legend('方位角误差', '俯仰角误差', '总误差', 'Location', 'northeast');
    grid on;
    
    % 理论分辨率分析
    subplot(2,2,2);
    azimuth_resolution = 0.886 * lambda / (N * d) * 180/pi;
    elevation_resolution = 0.886 * lambda / (M * d) * 180/pi;
    
    resolution_data = [azimuth_resolution, elevation_resolution];
    bar(resolution_data);
    set(gca, 'XTickLabel', {'方位角分辨率', '俯仰角分辨率'});
    ylabel('理论分辨率 (度)');
    title('阵列理论分辨率');
    grid on;
    
    % 计算复杂度分析
    subplot(2,2,3);
    n_angles = length(azimuth_range) * length(elevation_range);
    complexity_metrics = [n_angles, total_elements^3, K*log2(K)];
    bar(complexity_metrics);
    set(gca, 'XTickLabel', {'所提方法', 'MUSIC', 'FFT-based'});
    ylabel('计算复杂度指标');
    title('算法复杂度比较');
    grid on;
    
    % 性能总结
    subplot(2,2,4);
    text(0.1, 0.8, sprintf('阵列规模: %d×%d', M, N), 'FontSize', 12);
    text(0.1, 0.6, sprintf('最佳性能: %.2f°', min(results(:,4))), 'FontSize', 12);
    text(0.1, 0.4, sprintf('理论极限: %.2f°', min(azimuth_resolution, elevation_resolution)), 'FontSize', 12);
    text(0.1, 0.2, sprintf('接近度: %.1f%%', min(results(:,4))/min(azimuth_resolution, elevation_resolution)*100), 'FontSize', 12);
    axis off;
    title('性能总结');
    
    fprintf('性能分析完成\n');
end