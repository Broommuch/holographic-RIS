function result = verify_symbol_domain_gs(s_true, b_ref, A, t0)
%VERIFY_SYMBOL_DOMAIN_GS 验证 biased_gs_algorithm 在符号域是否能恢复QPSK符号
%
% 模型：
%   z = abs(A' * s_true + b_ref)
%
% 输入：
%   s_true : K x 1 真实QPSK符号向量
%   b_ref  : N x 1 已知参考信号，或标量参考信号
%   A      : K x N 已知信道/观测矩阵
%   t0     : GS算法迭代次数
%
% 输出：
%   result : 结构体，包含估计结果、硬判决结果、SER等信息
%
% 依赖：
%   biased_gs_algorithm(z, A, b_ref, t0)

    %% ================= 维度整理 =================
    s_true = s_true(:);

    [K, N] = size(A);

    if length(s_true) ~= K
        error('维度错误：s_true长度应等于size(A,1)。当前 length(s_true)=%d, size(A,1)=%d。', ...
              length(s_true), K);
    end

    % 如果 b_ref 是标量，则扩展成 N x 1
    if isscalar(b_ref)
        b_ref = b_ref * ones(N, 1);
    else
        b_ref = b_ref(:);
    end

    if length(b_ref) ~= N
        error('维度错误：b_ref长度应等于size(A,2)。当前 length(b_ref)=%d, size(A,2)=%d。', ...
              length(b_ref), N);
    end

    %% ================= 构造符号域接收幅度 =================
    % 复数叠加信号
    y_complex = A' * s_true + b_ref;

    % GS算法输入是幅度，不是功率
    z_obs = abs(y_complex);

    %% ================= 调用 biased GS 算法 =================
    s_est = biased_gs_algorithm(z_obs, A, b_ref, t0);

    s_est = s_est(:);

    %% ================= QPSK硬判决 =================
    qpsk_const = qammod((0:3).', 4, 'gray', 'UnitAveragePower', true);

    s_detect = zeros(K, 1);
    idx_detect = zeros(K, 1);

    for k = 1:K
        [~, idx_min] = min(abs(s_est(k) - qpsk_const));
        s_detect(k) = qpsk_const(idx_min);
        idx_detect(k) = idx_min - 1;
    end

    %% ================= 与真实QPSK符号对比 =================
    symbol_error_vec = (s_detect ~= s_true);
    num_error = sum(symbol_error_vec);
    SER = num_error / K;

    %% ================= 幅度重构误差 =================
    z_recon = abs(A' * s_est + b_ref);

    mag_mse = mean(abs(z_recon - z_obs).^2);
    mag_nmse = mag_mse / mean(abs(z_obs).^2);

    %% ================= 连续估计误差 =================
    est_mse = mean(abs(s_est - s_true).^2);

    %% ================= 输出结构体 =================
    result.s_true = s_true;
    result.b_ref = b_ref;
    result.A = A;
    result.z_obs = z_obs;
    result.s_est = s_est;
    result.s_detect = s_detect;
    result.idx_detect = idx_detect;
    result.symbol_error_vec = symbol_error_vec;
    result.num_error = num_error;
    result.SER = SER;
    result.z_recon = z_recon;
    result.mag_mse = mag_mse;
    result.mag_nmse = mag_nmse;
    result.est_mse = est_mse;

    %% ================= 打印结果 =================
    fprintf('\n========== Symbol-domain GS Verification ==========\n');
    fprintf('K = %d, N = %d\n', K, N);
    fprintf('GS iterations t0 = %d\n', t0);
    fprintf('Symbol errors = %d / %d\n', num_error, K);
    fprintf('SER = %.4f\n', SER);
    fprintf('Magnitude NMSE = %.4e\n', mag_nmse);
    fprintf('Continuous symbol MSE = %.4e\n', est_mse);

    if num_error == 0
        fprintf('Result: 所有QPSK符号均恢复正确。\n');
    else
        fprintf('Result: 存在符号恢复错误。\n');
    end

    %% ================= 可视化 =================
    figure;
    plot(real(s_true), imag(s_true), 'o', 'LineWidth', 1.5); hold on;
    plot(real(s_est), imag(s_est), 'x', 'LineWidth', 1.5);
    plot(real(s_detect), imag(s_detect), 's', 'LineWidth', 1.2);
    grid on;
    axis equal;
    xlabel('In-phase');
    ylabel('Quadrature');
    legend('True QPSK symbols', 'GS estimated symbols', 'Hard-decided symbols');
    title('Symbol-domain GS Recovery');

    figure;
    stem(1:N, z_obs, 'LineWidth', 1.4); hold on;
    stem(1:N, z_recon, '--', 'LineWidth', 1.2);
    grid on;
    xlabel('Observation index');
    ylabel('Magnitude');
    legend('Observed |A^H s + b|', 'Reconstructed |A^H \hat{s} + b|');
    title('Magnitude Fitting in Symbol Domain');

    figure;
    stem(1:K, symbol_error_vec, 'LineWidth', 1.4);
    grid on;
    xlabel('Symbol index');
    ylabel('Error flag');
    title('QPSK Symbol Error Positions');

end