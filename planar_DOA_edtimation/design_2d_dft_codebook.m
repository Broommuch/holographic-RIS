%% 3. 二维DFT码本设计
function codebook = design_2d_dft_codebook(M, N, K, bit_level)
% 设计二维DFT码本
% M, N: 阵列维度
% K: 码本数量
% bit_level: 量化比特数
    
    total_elements = M * N;
    codebook = zeros(total_elements, K);
    
    % 生成二维空间频率网格
    kx_values = linspace(-pi, pi, sqrt(K));
    ky_values = linspace(-pi, pi, sqrt(K));
    
    index = 1;
    for kx = kx_values(1:sqrt(K))
        for ky = ky_values(1:sqrt(K))
            % 生成二维DFT向量
            dft_vector = zeros(total_elements, 1);
            elem_index = 1;
            
            for m = 0:M-1
                for n = 0:N-1
                    phase = kx * n + ky * m;
                    dft_vector(elem_index) = exp(1j * phase);
                    elem_index = elem_index + 1;
                end
            end
            
            % 相位量化
            if bit_level < 32
                quantization_levels = 2^bit_level;
                phase_step = 2*pi / quantization_levels;
                
                for i = 1:total_elements
                    phase_continuous = angle(dft_vector(i));
                    quantized_phase = round(phase_continuous / phase_step) * phase_step;
                    dft_vector(i) = exp(1j * quantized_phase);
                end
            end
            
            codebook(:, index) = dft_vector / sqrt(total_elements);
            index = index + 1;
            
            if index > K
                break;
            end
        end
        if index > K
            break;
        end
    end
end