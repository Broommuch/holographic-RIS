%% 2. 不同比特数码本生成函数
function codebook = generate_quantized_codebook(N_ris, K, n_bits)
% 生成指定比特数的量化码本
% n_bits: 相位量化比特数
    
    % 生成基础DFT码本（连续相位）
    codebook_continuous = zeros(N_ris, K);
    for k = 0:K-1
        f_spatial = k * 2*pi / K;
%         f_spatial = (k-1 - floor(K/2)) * 2*pi / K;
        codebook_continuous(:, k+1) = exp(1j * f_spatial * (0:N_ris-1)') / sqrt(N_ris);
    end
    
    % 相位量化
    if n_bits >= 32 % 32比特以上视为连续相位
        codebook = codebook_continuous;
    else
        % 计算量化级别
        quantization_levels = 2^n_bits;
        phase_step = 2*pi / quantization_levels;
        
        % 量化相位
        codebook = zeros(N_ris, K);
        for k = 1:K
            for n = 1:N_ris
                phase_continuous = angle(codebook_continuous(n, k));
                % 量化到最近的离散相位
                quantized_phase = round(phase_continuous / phase_step) * phase_step;
                codebook(n, k) = exp(1j * quantized_phase);
            end
        end
        % 保持功率归一化
        codebook = codebook / sqrt(N_ris);
    end
end 