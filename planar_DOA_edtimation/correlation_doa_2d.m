%% 5. 二维相关性DOA估计算法
function [azimuth_est, elevation_est, correlation_map] = ...
    correlation_doa_2d(P_measured, codebook, azimuth_range, elevation_range, M, N, d, lambda)
% 二维DOA估计基于相关性分析
    
    n_azimuth = length(azimuth_range);
    n_elevation = length(elevation_range);
    correlation_map = zeros(n_elevation, n_azimuth);
    
    % 归一化测量功率
    P_measured_norm = P_measured / norm(P_measured);
    
    % 遍历所有角度组合
    for i = 1:n_elevation
        elevation = elevation_range(i);
        
        for j = 1:n_azimuth
            azimuth = azimuth_range(j);
            
            % 计算理论功率模式
            a_test = steering_vector_2d(azimuth, elevation, M, N, d, lambda);
            P_theoretical = zeros(size(P_measured));
            
            for k = 1:length(P_measured)
                w_k = codebook(:, k);
                P_theoretical(k) = abs(w_k' * a_test)^2;
            end
            
            % 归一化理论功率
            P_theoretical_norm = P_theoretical / norm(P_theoretical);
            
            % 计算相关性
            correlation_map(i, j) = P_measured_norm' * P_theoretical_norm;
        end
    end
    
    % 找到最大相关性位置
    [max_elev_idx, max_az_idx] = find(correlation_map == max(correlation_map(:)));
    azimuth_est = azimuth_range(max_az_idx(1));
    elevation_est = elevation_range(max_elev_idx(1));
end