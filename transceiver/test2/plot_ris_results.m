function plot_ris_results( ...
    SNRdB, ...
    MSE_theta, MSE_phi, ...
    crb_theta, crb_phi, ...
    SER, ...
    Ny, Nz, ...
    y_snapshot, ...
    fig_prefix)

if nargin < 10
    fig_prefix = '';
end

%% ================= Figure 1: DOA MSE vs SNR =================
figure;
semilogy(SNRdB, MSE_theta, '-o', 'LineWidth', 1.5); hold on;
semilogy(SNRdB, MSE_phi,   '-s', 'LineWidth', 1.5);
semilogy(SNRdB, crb_theta, '--', 'LineWidth', 1.5);
semilogy(SNRdB, crb_phi,   '-.', 'LineWidth', 1.5);
grid on;

xlabel('SNR (dB)');
ylabel('Mean Squared Error');
legend( ...
    'Elevation MSE (GN)', ...
    'Azimuth MSE (GN)', ...
    'Elevation CRB', ...
    'Azimuth CRB', ...
    'Location', 'southwest');

title('DOA Estimation Performance');

if ~isempty(fig_prefix)
    saveas(gcf, [fig_prefix '_DOA_MSE.png']);
end


%% ================= Figure 2: SER vs SNR =================
figure;
semilogy(SNRdB, SER, '-o', 'LineWidth', 1.5);
grid on;

xlabel('SNR (dB)');
ylabel('Symbol Error Rate');
title('QPSK Symbol Error Rate');

if ~isempty(fig_prefix)
    saveas(gcf, [fig_prefix '_SER.png']);
end


%% ================= Figure 3: RIS Power Heatmap =================
% y_snapshot is assumed to be (V*A*s + b) or noisy version

% %功率图绘制默认
% if ~isempty(y_snapshot)
% 
%     P = abs(y_snapshot).^2;
%     Pmap = reshape(P, Nz, Ny);
% 
%     figure;
%     imagesc(Pmap);
%     axis equal tight;
%     colormap jet;
%     colorbar;
% 
%     xlabel('RIS element index (y-axis)');
%     ylabel('RIS element index (z-axis)');
%     title('RIS Power Distribution');
% 
%     if ~isempty(fig_prefix)
%         saveas(gcf, [fig_prefix '_RIS_PowerMap.png']);
%     end
% end

% %功率图绘制图像插值
% if ~isempty(y_snapshot)
%     P = abs(y_snapshot).^2;
%     Pmap = reshape(P, Nz, Ny);
%     
%     % 插值平滑处理
%     [Ny_grid, Nz_grid] = size(Pmap);
%     [X, Y] = meshgrid(1:Ny, 1:Nz);
%     [Xi, Yi] = meshgrid(linspace(1, Ny, 5*Ny), linspace(1, Nz, 5*Nz));
%     Pmap_interp = interp2(X, Y, Pmap, Xi, Yi, 'spline');
%     
%     figure;
%     imagesc(Pmap_interp);
%     axis equal tight;
%     colormap(jet);
%     colorbar;
%     
%     xlabel('RIS element index (y-axis)');
%     ylabel('RIS element index (z-axis)');
%     title('RIS Power Distribution (Interpolated)');
%     
%     if ~isempty(fig_prefix)
%         saveas(gcf, [fig_prefix '_RIS_PowerMap_Smooth.png']);
%     end
% end

%功率图绘制 - 高密度插值 + 3D可视化
if ~isempty(y_snapshot)
    P = abs(y_snapshot).^2;
    Pmap = reshape(P, Nz, Ny);
    
    % 1. 增大插值密度（可根据需要调整插值倍数）
    interp_factor = 4;  % 插值倍数，可根据需要调整，如8, 10, 15等
    
    [Nz_orig, Ny_orig] = size(Pmap);
    [X, Y] = meshgrid(1:Ny_orig, 1:Nz_orig);
    
    % 创建高密度插值网格
    Ny_interp = interp_factor * Ny_orig;
    Nz_interp = interp_factor * Nz_orig;
    [Xi, Yi] = meshgrid(linspace(1, Ny_orig, Ny_interp), ...
                        linspace(1, Nz_orig, Nz_interp));
    
    % 进行插值（可以尝试不同的插值方法）
    Pmap_interp = interp2(X, Y, Pmap, Xi, Yi, 'cubic');  % 用cubic插值更平滑
    
    % 2. 2D高密度插值图
    figure('Position', [100, 100, 1200, 500]);
    
    subplot(1, 2, 1);
    imagesc(Pmap_interp);
    axis equal tight;
    colormap(jet);
    colorbar;
    title(['2D RIS Power Distribution (Interp Factor: ', num2str(interp_factor), 'x)']);
    xlabel('RIS element index (y-axis)');
    ylabel('RIS element index (z-axis)');
    
%     % 添加原始数据点标记
%     hold on;
%     scatter(X(:), Y(:), 30, 'k', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1);
%     hold off;
    
    % 3. 3D可视化
    subplot(1, 2, 2);
    
    % 创建3D网格
    [X3d, Y3d] = meshgrid(1:size(Pmap_interp, 2), 1:size(Pmap_interp, 1));
    
    % 绘制3D曲面
    surf(X3d, Y3d, Pmap_interp, 'EdgeColor', 'none', 'FaceColor', 'interp');
    
    % 设置视角
    view(45, 30);  % 方位角45°，仰角30°
    
    % 设置颜色和光照
    colormap(jet);
    shading interp;
    
    % 添加光照效果
    light('Position', [1, 1, 1], 'Style', 'infinite');
    lighting gouraud;
    material([0.3, 0.8, 0.3, 10, 0.5]);
    
    axis tight;
    grid on;
    colorbar;
    
    xlabel('RIS element index (y-axis)');
    ylabel('RIS element index (z-axis)');
    zlabel('Power');
    title('3D RIS Power Distribution');
    
    % 添加颜色图例
    cb = colorbar;
    ylabel(cb, 'Power (dB)');
    
    if ~isempty(fig_prefix)
        saveas(gcf, [fig_prefix '_RIS_PowerMap_Interp3D.png']);
    end
end

% %功率图绘制pcolor
% if ~isempty(y_snapshot)
%     P = abs(y_snapshot).^2;
%     Pmap = reshape(P, Nz, Ny);
%     
%     figure;
%     % 创建扩展的网格
%     [Ny_grid, Nz_grid] = size(Pmap);
%     [X, Y] = meshgrid(1:Ny_grid+1, 1:Nz_grid+1);
%     
%     % 扩展数据矩阵
%     Pmap_expanded = zeros(Nz_grid+1, Ny_grid+1);
%     Pmap_expanded(1:Nz_grid, 1:Ny_grid) = Pmap;
%     
%     pcolor(X, Y, Pmap_expanded);
%     shading interp;  % 关键：插值着色
%     axis equal tight;
%     colormap(jet);
%     colorbar;
%     
%     xlabel('RIS element index (y-axis)');
%     ylabel('RIS element index (z-axis)');
%     title('RIS Power Distribution');
%     
%     if ~isempty(fig_prefix)
%         saveas(gcf, [fig_prefix '_RIS_PowerMap_Smooth.png']);
%     end
% end

% %高斯滤波平滑
% if ~isempty(y_snapshot)
%     P = abs(y_snapshot).^2;
%     Pmap = reshape(P, Nz, Ny);
%     
%     % 高斯滤波平滑
%     sigma = 0.8;  % 平滑系数
%     Pmap_smooth = imgaussfilt(Pmap, sigma);
%     
%     figure;
%     imagesc(Pmap_smooth);
%     axis equal tight;
%     colormap(jet);
%     colorbar;
%     
%     xlabel('RIS element index (y-axis)');
%     ylabel('RIS element index (z-axis)');
%     title('RIS Power Distribution (Smoothed)');
%     
%     if ~isempty(fig_prefix)
%         saveas(gcf, [fig_prefix '_RIS_PowerMap_Smooth.png']);
%     end
% end

% %surf函数
% if ~isempty(y_snapshot)
%     P = abs(y_snapshot).^2;
%     Pmap = reshape(P, Nz, Ny);
%     
%     figure;
%     surf(Pmap, 'EdgeColor', 'none', 'FaceColor', 'interp');
%     view(2);  % 俯视图
%     axis equal tight;
%     colormap(jet);
%     colorbar;
%     
%     xlabel('RIS element index (y-axis)');
%     ylabel('RIS element index (z-axis)');
%     title('RIS Power Distribution');
%     
%     if ~isempty(fig_prefix)
%         saveas(gcf, [fig_prefix '_RIS_PowerMap_Smooth.png']);
%     end
% end

% % 等高线图
% if ~isempty(y_snapshot)
%     P = abs(y_snapshot).^2;
%     Pmap = reshape(P, Nz, Ny);
%     
%     figure;
%     % 等高线填充
%     contourf(Pmap, 20, 'LineStyle', 'none');
%     axis equal tight;
%     colormap(jet);
%     colorbar;
%     
%     xlabel('RIS element index (y-axis)');
%     ylabel('RIS element index (z-axis)');
%     title('RIS Power Distribution (Contour)');
%     
%     if ~isempty(fig_prefix)
%         saveas(gcf, [fig_prefix '_RIS_PowerMap_Contour.png']);
%     end
% end

% % 3D视图
% if ~isempty(y_snapshot)
%     P = abs(y_snapshot).^2;
%     Pmap = reshape(P, Nz, Ny);
%     
%     figure;
%     surf(Pmap, 'EdgeColor', 'none', 'FaceColor', 'interp');
%     % 删除或注释掉 view(2) 以显示3D视图
%     % view(2);  % 注释掉这行
%     axis tight;
%     colormap(jet);
%     colorbar;
%     
%     xlabel('RIS element index (y-axis)');
%     ylabel('RIS element index (z-axis)');
%     zlabel('Power');
%     title('3D RIS Power Distribution');
%     
%     % 添加光照和材质效果增强3D感
%     light;
%     lighting gouraud;
%     material dull;
%     
%     if ~isempty(fig_prefix)
%         saveas(gcf, [fig_prefix '_RIS_PowerMap_3D.png']);
%     end
% end


end
