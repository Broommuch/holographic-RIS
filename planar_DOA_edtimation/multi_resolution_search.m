%% 9. 高级功能：多分辨率搜索
function [az_refined, el_refined] = multi_resolution_search(...
    P_measured, codebook, az_coarse, el_coarse, M, N, d, lambda)
% 多分辨率角度搜索
    
    % 第一轮：粗搜索
    az_range_coarse = az_coarse-10:5:az_coarse+10;
    el_range_coarse = max(0, el_coarse-10):5:min(60, el_coarse+10);
    
    [az_mid, el_mid, ~] = correlation_doa_2d(P_measured, codebook, ...
                                            az_range_coarse, el_range_coarse, M, N, d, lambda);
    
    % 第二轮：精搜索
    az_range_fine = az_mid-2:0.5:az_mid+2;
    el_range_fine = max(0, el_mid-2):0.5:min(60, el_mid+2);
    
    [az_refined, el_refined, ~] = correlation_doa_2d(P_measured, codebook, ...
                                                    az_range_fine, el_range_fine, M, N, d, lambda);
end