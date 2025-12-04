function codebook = design_optimal_codebook(N_ris, K, design_type)
% 设计最优空间频率编码码本
% 输入: N_ris - RIS单元数, K - 码本数量, design_type - 码本类型
% 输出: codebook - N_ris × K 的码本矩阵

    switch design_type
        case 'dft'  % 离散傅里叶变换码本（最优选择）
            codebook = zeros(N_ris, K);
            for k = 1:K
                % 生成DFT向量，覆盖完整空间频率范围
                f_spatial = (k-1 - floor(K/2)) * 2*pi / K;
                codebook(:, k) = exp(1j * f_spatial * (0:N_ris-1)') / sqrt(N_ris);
            end
            
        case 'hadamard'  % 哈达玛码本（正交性好）
            H = hadamard(K);
            % 取前N_ris行，映射到相位
            codebook = H(1:N_ris, 1:K);
            codebook = exp(1j * pi * (1 - codebook)/2);
            
        case 'random_structured'  % 结构化随机码本
            rng(42); % 可重复性
            base_phases = linspace(0, 2*pi, K);
            codebook = zeros(N_ris, K);
            for k = 1:K
                % 在基础相位上添加小随机扰动
                phases = base_phases(k) + 0.1 * randn(N_ris, 1);
                codebook(:, k) = exp(1j * phases);
            end
            
        case 'sparse'  % 稀疏码本（压缩感知优化）
            codebook = zeros(N_ris, K);
            for k = 1:K
                % 每个码本只激活部分单元
                active_ratio = 0.3; % 30%单元激活
                active_mask = rand(N_ris, 1) < active_ratio;
                phases = 2*pi * rand(N_ris, 1);
                codebook(:, k) = active_mask .* exp(1j * phases);
            end
    end
end