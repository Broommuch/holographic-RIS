% 算法1: 基于FFT的空间频率估计
function theta_est = fft_based_doa(P_sequence, codebook, angles, lambda, d)
    K = length(P_sequence);
    N_ris = size(codebook, 1);
    
    % FFT分析
    P_fft = fft(P_sequence - mean(P_sequence));
    P_fft_mag = abs(P_fft(1:floor(K/2)+1));
    
    % 空间频率映射
    f_spatial_axis = (0:length(P_fft_mag)-1) / K;
    theta_fft_axis = asind(f_spatial_axis * lambda / d);
    
    % 峰值检测
    [~, peak_idx] = max(P_fft_mag(2:end)); % 跳过直流
    theta_est = theta_fft_axis(peak_idx + 1);
end