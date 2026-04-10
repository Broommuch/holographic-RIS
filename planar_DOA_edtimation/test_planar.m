%% 二维平面RIS阵列DOA估计系统
% 基于DFT码本和相关性算法的二维角度估计
clear; clc; close all;

%% 1. 系统参数设置
fprintf('=== 二维平面RIS阵列DOA估计系统 ===\n');

% 物理参数
f0 = 3.5e9;        % 载波频率 3.5GHz
c = 3e8;          % 光速
lambda = c/f0;    % 波长
d = lambda/2;     % 单元间距

% 二维阵列参数
M = 8;            % 行数 (y方向)
N = 8;            % 列数 (x方向)
total_elements = M * N;

% 码本参数
K = 16;           % 码本数量
SNR_dB = 20;      % 信噪比

% 角度搜索范围
azimuth_range = -60:2:60;    % 方位角搜索范围 (度)
elevation_range = 0:2:60;    % 俯仰角搜索范围 (度)

% 真实入射角度
azimuth_true = 25;    % 方位角 (度)
elevation_true = 30;  % 俯仰角 (度)

fprintf('阵列规模: %d×%d = %d 单元\n', M, N, total_elements);
fprintf('真实入射方向: 方位角%.1f°, 俯仰角%.1f°\n', azimuth_true, elevation_true);




%% 6. 主程序执行
fprintf('\n=== 开始二维DOA估计 ===\n');

% 生成码本
codebook_2d = design_2d_dft_codebook(M, N, K, 4); % 使用4-bit量化

% 模拟功率测量
P_measured_2d = simulate_power_measurement_2d(codebook_2d, azimuth_true, ...
                                            elevation_true, M, N, d, lambda, SNR_dB);

% 执行二维DOA估计
tic;
[azimuth_est, elevation_est, correlation_map] = ...
    correlation_doa_2d(P_measured_2d, codebook_2d, azimuth_range, ...
                      elevation_range, M, N, d, lambda);
estimation_time = toc;

% 计算估计误差
azimuth_error = abs(azimuth_est - azimuth_true);
elevation_error = abs(elevation_est - elevation_true);
total_error = sqrt(azimuth_error^2 + elevation_error^2);

fprintf('估计结果:\n');
fprintf('  方位角: 真实%.1f° -> 估计%.1f° (误差: %.2f°)\n', ...
        azimuth_true, azimuth_est, azimuth_error);
fprintf('  俯仰角: 真实%.1f° -> 估计%.1f° (误差: %.2f°)\n', ...
        elevation_true, elevation_est, elevation_error);
fprintf('  总误差: %.2f°\n', total_error);
fprintf('  计算时间: %.3f秒\n', estimation_time);

%% 7. 结果可视化
fprintf('\n=== 生成可视化结果 ===\n');

% 图1: 相关性图谱
figure('Position', [100, 100, 1200, 800]);
subplot(2,3,1);
imagesc(azimuth_range, elevation_range, correlation_map);
xlabel('方位角 (度)');
ylabel('俯仰角 (度)');
title('二维相关性图谱');
colorbar;
hold on;
plot(azimuth_true, elevation_true, 'rx', 'MarkerSize', 15, 'LineWidth', 3);
plot(azimuth_est, elevation_est, 'go', 'MarkerSize', 10, 'LineWidth', 2);
legend('真实角度', '估计角度', 'Location', 'best');

% 图2: 三维相关性曲面
subplot(2,3,2);
[Az, El] = meshgrid(azimuth_range, elevation_range);
surf(Az, El, correlation_map, 'EdgeColor', 'none');
xlabel('方位角 (度)');
ylabel('俯仰角 (度)');
zlabel('相关性');
title('三维相关性曲面');
colorbar;

% 图3: 功率测量序列
subplot(2,3,3);
plot(1:K, P_measured_2d, 'b-o', 'LineWidth', 1.5);
xlabel('码本索引');
ylabel('测量功率');
title('功率测量序列');
grid on;

% 图4: 阵列模式示例
subplot(2,3,4);
% 计算阵列方向图
az_test = -90:2:90;
el_test = 0:2:90;
pattern = zeros(length(el_test), length(az_test));

for i = 1:length(el_test)
    for j = 1:length(az_test)
        a_temp = steering_vector_2d(az_test(j), el_test(i), M, N, d, lambda);
        % 使用第一个码本
        pattern(i, j) = abs(codebook_2d(:,1)' * a_temp)^2;
    end
end

imagesc(az_test, el_test, 10*log10(pattern/max(pattern(:))));
xlabel('方位角 (度)');
ylabel('俯仰角 (度)');
title('阵列方向图 (dB)');
colorbar;

% 图5: 误差分析
subplot(2,3,5);
errors = [azimuth_error, elevation_error, total_error];
bar(errors);
set(gca, 'XTickLabel', {'方位角误差', '俯仰角误差', '总误差'});
ylabel('误差 (度)');
title('估计误差分析');
grid on;

% 图6: 码本相位分布
subplot(2,3,6);
% 显示第一个码本的相位分布
codebook_phase = angle(codebook_2d(:,1));
codebook_phase_2d = reshape(codebook_phase, M, N);
imagesc(codebook_phase_2d);
colorbar;
xlabel('列索引');
ylabel('行索引');
title('码本相位分布 (弧度)');



% 执行性能分析
analyze_performance_2d(codebook_2d,M,N,d,lambda,azimuth_range,...
    elevation_range,total_elements,K);



% 测试多分辨率搜索
fprintf('\n=== 多分辨率搜索测试 ===\n');
[az_refined, el_refined] = multi_resolution_search(...
    P_measured_2d, codebook_2d, azimuth_est, elevation_est, M, N, d, lambda);

fprintf('粗估计: 方位角%.1f°, 俯仰角%.1f°\n', azimuth_est, elevation_est);
fprintf('精估计: 方位角%.1f°, 俯仰角%.1f°\n', az_refined, el_refined);
fprintf('精度提升: %.2f°\n', total_error - sqrt((az_refined-azimuth_true)^2 + (el_refined-elevation_true)^2));

%% 10. 系统验证与保存结果
fprintf('\n=== 系统验证完成 ===\n');

% 保存重要结果
save('2D_RIS_DOA_Results.mat', 'azimuth_est', 'elevation_est', ...
     'correlation_map', 'codebook_2d', 'P_measured_2d');

% 生成性能报告
report_file = '2D_DOA_Performance_Report.txt';
fid = fopen(report_file, 'w');
fprintf(fid, '二维RIS阵列DOA估计性能报告\n');
fprintf(fid, '生成时间: %s\n\n', datestr(now));

fprintf(fid, '系统参数:\n');
fprintf(fid, '阵列规模: %d×%d\n', M, N);
fprintf(fid, '码本数量: %d\n', K);
fprintf(fid, '信噪比: %d dB\n\n', SNR_dB);

fprintf(fid, '估计结果:\n');
fprintf(fid, '真实角度: 方位角%.1f°, 俯仰角%.1f°\n', azimuth_true, elevation_true);
fprintf(fid, '估计角度: 方位角%.1f°, 俯仰角%.1f°\n', azimuth_est, elevation_est);
fprintf(fid, '估计误差: %.2f°\n\n', total_error);

fprintf(fid, '性能指标:\n');
fprintf(fid, '计算时间: %.3f秒\n', estimation_time);
fprintf(fid, '理论分辨率: 方位角%.2f°, 俯仰角%.2f°\n', ...
        0.886*lambda/(N*d)*180/pi, 0.886*lambda/(M*d)*180/pi);

fclose(fid);
fprintf('性能报告已保存: %s\n', report_file);

fprintf('\n=== 二维DOA估计系统运行完成 ===\n');
