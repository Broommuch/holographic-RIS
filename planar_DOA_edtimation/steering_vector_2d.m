%% 2. 二维导向矢量函数
function a = steering_vector_2d(azimuth, elevation, M, N, d, lambda)
% 计算二维平面阵列的导向矢量
% azimuth: 方位角 (度)
% elevation: 俯仰角 (度)
% M: 行数, N: 列数
    
    azimuth_rad = deg2rad(azimuth);
    elevation_rad = deg2rad(elevation);
    
    % 计算波达方向向量
    k_vector = 2*pi/lambda * [sin(elevation_rad)*cos(azimuth_rad);
                             sin(elevation_rad)*sin(azimuth_rad);
                             cos(elevation_rad)];
    
    a = zeros(M*N, 1);
    index = 1;
    
    for m = 0:M-1
        for n = 0:N-1
            % 单元位置 (x, y, z)
            r = [n*d; m*d; 0];
            
            % 相位延迟
            phase_delay = k_vector' * r;
            a(index) = exp(1j * phase_delay);
            index = index + 1;
        end
    end
end