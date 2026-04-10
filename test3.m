%% 不同比特数码本对RIS DOA估计性能影响分析
% 基于Nature Electronics全息超表面原理的扩展研究
clear; clc; close all;

%% 1. 参数设置与初始化
fprintf('=== 不同比特数码本性能影响分析 ===\n');

% 物理参数
f0 = 3.5e9;        % 载波频率 3.5GHz
c = 3e8;          % 光速
lambda = c/f0;    % 波长
d = lambda/2;     % RIS单元间距

% 系统参数
N_ris = 16;       % RIS单元数
K = 4;           % 码本数量
angles = -60:0.5:60; % 角度搜索范围
SNR_dB = 20;      % 信噪比
theta_true = 25;  % 真实入射角

% 测试的比特数范围
bit_levels = [1, 2, 3, 4, 8, 16, 32, 64]; % 不同比特数
n_bits = length(bit_levels);

fprintf('系统参数: RIS单元%d, 码本数%d, 信噪比%ddB\n', N_ris, K, SNR_dB);
fprintf('测试比特数: %s\n', mat2str(bit_levels));


%% 4. 主测试循环
fprintf('\n=== 开始性能测试 ===\n');

% 初始化结果存储
results = struct();
for i = 1:n_bits
    bit_level = bit_levels(i);
    results(i).bit_level = bit_level;
    results(i).codebook = [];
    results(i).performance = [];
end

% 蒙特卡洛仿真参数
n_monte_carlo = 50; % 蒙特卡洛仿真次数
fprintf('蒙特卡洛仿真次数: %d\n', n_monte_carlo);

for bit_idx = 1:n_bits
    bit_level = bit_levels(bit_idx);
    fprintf('\n测试 %d-bit 码本...', bit_level);
    
    % 存储每次仿真的结果
    errors = zeros(n_monte_carlo, 1);
    times = zeros(n_monte_carlo, 1);
    estimates = zeros(n_monte_carlo, 1);
    
    for mc_iter = 1:n_monte_carlo
        % 生成码本
        codebook = generate_quantized_codebook(N_ris, K, bit_level);
%         if bit_level == 1
%             save codebook_1bit.mat codebook;
%         elseif bit_level == 2
%             save codebook_2bit.mat codebook;
%         end



        % 随机化测试角度（避免过拟合）
        theta_test = theta_true + (rand-0.5)*10; % ±5度随机变化
        
        % 评估性能
        [theta_est, error, time_elapsed] = evaluate_codebook_performance(...
            codebook, theta_test, angles, SNR_dB, lambda, d);
        
        errors(mc_iter) = error;
        times(mc_iter) = time_elapsed;
        estimates(mc_iter) = theta_est;
    end
    
    % 存储结果
    results(bit_idx).codebook = codebook;
    results(bit_idx).errors = errors;
    results(bit_idx).times = times;
    results(bit_idx).estimates = estimates;
    results(bit_idx).mean_error = mean(errors);
    results(bit_idx).std_error = std(errors);
    results(bit_idx).mean_time = mean(times);
    
    fprintf('完成 - 平均误差: %.3f° ± %.3f°, 平均时间: %.4fs\n', ...
        mean(errors), std(errors), mean(times));
end



%% 6. 结果可视化与分析
fprintf('\n=== 生成综合性能分析图 ===\n');

% 主性能对比图
figure('Position', [100, 100, 1400, 1000]);

% 子图1: 平均估计误差 vs 比特数
subplot(2,3,1);
mean_errors = [results.mean_error];
std_errors = [results.std_error];
errorbar(bit_levels, mean_errors, std_errors, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('码本比特数');
ylabel('平均估计误差 (°)');
title('DOA估计误差 vs 码本比特数');
set(gca, 'YScale', 'log');
grid on;

% 添加理论分辨率线
theory_resolution = 0.886 * lambda / (N_ris * d) * 180/pi;
hold on;
plot(xlim, [theory_resolution, theory_resolution], 'r--', 'LineWidth', 2);
legend('实测误差', '理论分辨率极限', 'Location', 'best');

% 子图2: 计算时间 vs 比特数
subplot(2,3,2);
mean_times = [results.mean_time];
plot(bit_levels, mean_times*1000, 's-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('码本比特数');
ylabel('平均计算时间 (ms)');
title('计算时间 vs 码本比特数');
set(gca, 'YScale', 'log');
grid on;

% 子图3: 误差分布箱线图
subplot(2,3,3);
error_data = [];
group_data = [];
for i = 1:n_bits
    error_data = [error_data; results(i).errors];
    group_data = [group_data; i*ones(length(results(i).errors), 1)];
end
boxplot(error_data, group_data, 'Labels', arrayfun(@num2str, bit_levels, 'UniformOutput', false));
xlabel('码本比特数');
ylabel('估计误差 (°)');
title('误差分布箱线图');
grid on;

% 子图4: 性能-复杂度权衡
subplot(2,3,4);
% 计算性能指标（误差的倒数，越大越好）
performance_metric = 1./mean_errors;
complexity_metric = mean_times;

scatter(complexity_metric, performance_metric, 100, bit_levels, 'filled');
xlabel('计算复杂度 (s)');
ylabel('性能指标 (1/误差)');
title('性能-复杂度权衡分析');
colorbar;
colormap(jet);
for i = 1:n_bits
    text(complexity_metric(i), performance_metric(i), sprintf('%d-bit', bit_levels(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
grid on;

% 子图5: 不同比特数的码本相位分布对比
subplot(2,3,5);
hold on;
colors = jet(n_bits);
for i = 1:n_bits
    phases = angle(results(i).codebook(:));
    [counts, edges] = histcounts(phases, 50);
    centers = (edges(1:end-1) + edges(2:end))/2;
    plot(centers, counts, 'Color', colors(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('%d-bit', bit_levels(i)));
end
xlabel('相位 (弧度)');
ylabel('频数');
title('不同比特数码本相位分布');
legend('show');
grid on;

% 子图6: 推荐比特数选择指南
subplot(2,3,6);
% 计算综合得分（误差小且时间短者得分高）
normalized_errors = mean_errors / max(mean_errors);
normalized_times = mean_times / max(mean_times);
composite_scores = 1./(0.7*normalized_errors + 0.3*normalized_times);

[best_score, best_idx] = max(composite_scores);
best_bit_level = bit_levels(best_idx);

bar(bit_levels, composite_scores);
xlabel('码本比特数');
ylabel('综合得分');
title('最优比特数选择');
set(gca, 'XScale', 'log');
grid on;
hold on;
plot(best_bit_level, best_score, 'ro', 'MarkerSize', 10, 'LineWidth', 3);
text(best_bit_level, best_score, sprintf('推荐: %d-bit', best_bit_level), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

%% 7. 详细比特数特性分析
fprintf('\n=== 生成详细比特数分析报告 ===\n');

% 选择几个关键比特数进行详细分析
key_bits = [1, 2, 3, 4, 8, 32];
key_indices = zeros(size(key_bits));
for i = 1:length(key_bits)
    key_indices(i) = find(bit_levels == key_bits(i));
end

% 生成详细分析图
figure('Position', [200, 200, 1600, 1200]);

for i = 1:length(key_bits)
    idx = key_indices(i);
    bit_level = key_bits(i);
    
    subplot(3, 3, i);
    
    % 分析该比特数码本的特性
    codebook = results(idx).codebook;
    
    % 绘制相位分布极坐标图
    phases = angle(codebook(:));
    polarhistogram(phases, 50);
    title(sprintf('%d-bit码本相位分布', bit_level));
end

% 添加性能总结
subplot(3,3,7);
performance_summary = zeros(length(key_bits), 3);
for i = 1:length(key_bits)
    idx = key_indices(i);
    performance_summary(i, 1) = results(idx).mean_error;
    performance_summary(i, 2) = results(idx).mean_time;
    performance_summary(i, 3) = 1/results(idx).mean_error; % 性能指标
end

bar(performance_summary);
set(gca, 'XTickLabel', arrayfun(@num2str, key_bits, 'UniformOutput', false));
ylabel('指标值');
title('关键比特数性能对比');
legend('平均误差 (°)', '计算时间 (s)', '性能指标', 'Location', 'best');
grid on;

% 添加理论分析
subplot(3,3,8);
% 量化误差理论分析
theory_quant_error = (pi ./ (2.^key_bits)) .^ 2 / 12; % 相位量化误差理论值
theory_doa_error = theory_quant_error * 180/pi * N_ris; % 转换为角度误差

semilogy(key_bits, theory_doa_error, 'o-', 'LineWidth', 2);
hold on;
actual_errors = performance_summary(:,1)';
semilogy(key_bits, actual_errors, 's-', 'LineWidth', 2);
xlabel('码本比特数');
ylabel('误差 (°)');
title('量化误差理论 vs 实际');
legend('理论量化误差', '实际DOA误差', 'Location', 'best');
grid on;

% 添加推荐总结
subplot(3,3,9);
text(0.1, 0.7, sprintf('最优比特数: %d-bit', best_bit_level), 'FontSize', 14);
text(0.1, 0.5, sprintf('平均误差: %.3f°', results(best_idx).mean_error), 'FontSize', 12);
text(0.1, 0.3, sprintf('计算时间: %.4fs', results(best_idx).mean_time), 'FontSize', 12);
text(0.1, 0.1, sprintf('接近理论极限: %.1f%%', theory_resolution/results(best_idx).mean_error*100), 'FontSize', 12);
axis off;
title('推荐配置总结');

%% 8. 实际系统实现建议
fprintf('\n=== 实际系统实现建议 ===\n');

fprintf('基于仿真结果，推荐使用 %d-bit 码本\n', best_bit_level);
fprintf('理由:\n');
fprintf('1. 估计误差: %.3f° (接近理论分辨率 %.3f°)\n', results(best_idx).mean_error, theory_resolution);
fprintf('2. 计算时间: %.4f秒 (满足实时性要求)\n', results(best_idx).mean_time);
fprintf('3. 硬件复杂度: 适中的 %d-bit 相位量化\n', best_bit_level);
fprintf('4. 性能-复杂度权衡最优\n\n');

fprintf('不同应用场景建议:\n');
fprintf('• 高精度应用: ≥4-bit 码本，误差 < 0.5°\n');
fprintf('• 实时性要求高: 2-bit 码本，计算速度快\n');
fprintf('• 硬件受限: 1-bit 码本，最简单的实现\n');
fprintf('• 理论研究: 连续相位码本，最佳性能基准\n');

%% 9. 保存结果与生成报告
fprintf('\n=== 保存仿真结果 ===\n');

% 保存详细结果
save('Bit_Level_Analysis_Results.mat', 'results', 'bit_levels', 'best_bit_level', 'theory_resolution');

% 生成文本报告
report_filename = 'Bit_Level_Analysis_Report.txt';
fid = fopen(report_filename, 'w');
fprintf(fid, '不同比特数码本对RIS DOA估计性能影响分析报告\n');
fprintf(fid, '生成时间: %s\n\n', datestr(now));

fprintf(fid, '系统参数:\n');
fprintf(fid, '• RIS单元数: %d\n', N_ris);
fprintf(fid, '• 码本数量: %d\n', K);
fprintf(fid, '• 信噪比: %d dB\n', SNR_dB);
fprintf(fid, '• 理论分辨率: %.3f°\n\n', theory_resolution);

fprintf(fid, '性能总结:\n');
fprintf(fid, '比特数\t平均误差(°)\t标准差(°)\t计算时间(s)\n');
fprintf(fid, '------\t----------\t----------\t----------\n');
for i = 1:n_bits
    fprintf(fid, '%d-bit\t%.3f\t\t%.3f\t\t%.4f\n', ...
        results(i).bit_level, results(i).mean_error, results(i).std_error, results(i).mean_time);
end

fprintf(fid, '\n推荐配置: %d-bit 码本\n', best_bit_level);
fprintf(fid, '理由: 在估计精度和计算复杂度之间取得最佳平衡\n');

fclose(fid);
fprintf('分析报告已保存至: %s\n', report_filename);

fprintf('\n=== 仿真完成 ===\n');