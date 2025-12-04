%% 5. 码本特性分析函数
function analyze_codebook_properties(codebook, bit_level)
% 分析码本的各种特性
    
    [N_ris, K] = size(codebook);
    
    % 1. 相位分布分析
    phases = angle(codebook(:));
    
    figure('Position', [100, 100, 1200, 800]);
    
    % 子图1: 相位分布直方图
    subplot(2,3,1);
    histogram(phases, 50);
    xlabel('相位 (弧度)');
    ylabel('频数');
    title(sprintf('%d-bit码本相位分布', bit_level));
    grid on;
    
    % 子图2: 自相关矩阵
    subplot(2,3,2);
    correlation_matrix = abs(codebook' * codebook);
    imagesc(correlation_matrix);
    colorbar;
    xlabel('码本索引');
    ylabel('码本索引');
    title('码本自相关矩阵');
    
    % 子图3: 正交性度量
    subplot(2,3,3);
    ideal_orthogonal = eye(K);
    orthogonality_error = norm(correlation_matrix - ideal_orthogonal, 'fro');
    bar(1, orthogonality_error);
    ylabel('正交性误差');
    title(sprintf('正交性误差: %.4f', orthogonality_error));
    grid on;
    
    % 子图4: 码本向量在复平面上的分布
    subplot(2,3,4);
    plot(real(codebook(:)), imag(codebook(:)), '.', 'MarkerSize', 1);
    xlabel('实部');
    ylabel('虚部');
    title('码本向量分布');
    axis equal;
    grid on;
    
    % 子图5: 空间频率响应
    subplot(2,3,5);
    angles_test = -60:1:60;
    response = zeros(length(angles_test), 1);
    steering_vector = @(theta) exp(1j * 2*pi * 0.5 * (0:N_ris-1)' * sind(theta));
    
    for i = 1:length(angles_test)
        a_theta = steering_vector(angles_test(i));
        % 使用第一个码本作为示例
        response(i) = abs(codebook(:,1)' * a_theta)^2;
    end
    
    plot(angles_test, 10*log10(response/max(response)));
    xlabel('角度 (°)');
    ylabel('归一化响应 (dB)');
    title('示例码本波束模式');
    grid on;
    
    % 子图6: 码本差异性
    subplot(2,3,6);
    diversity_metric = zeros(K-1, 1);
    for k = 1:K-1
        diversity_metric(k) = norm(codebook(:,k) - codebook(:,k+1));
    end
    plot(1:K-1, diversity_metric, 'o-');
    xlabel('相邻码本对');
    ylabel('差异度量');
    title('码本序列差异性');
    grid on;
end