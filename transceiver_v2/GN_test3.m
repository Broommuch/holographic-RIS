%% 这个函数测试了收敛速度的效果，说明对初值真的很敏感
clc; clear; close all;

%% ================= 1. 参数 =================
N = 128; d = 0.5;
theta_true = 30; phi_true = 20;
u_true = sin(deg2rad(theta_true))*cos(deg2rad(phi_true));
v_true = sin(deg2rad(theta_true))*sin(deg2rad(phi_true));
[x_idx, y_idx] = meshgrid(0:N-1, 0:N-1);
cx = (N-1)/2; cy = (N-1)/2;

%% ================= 2. 模拟测量 =================
lambda_wave = 633e-9; pixel_size = 5e-6;
R = exp(1i*2*pi/lambda_wave * ((x_idx-cx)*pixel_size*sin(deg2rad(20)) + (y_idx-cy)*pixel_size*sin(deg2rad(10))));
O_true = exp(1i * 2*pi*d*((x_idx-cx)*u_true + (y_idx-cy)*v_true));
I = abs(O_true + R).^2;

%% ================= 3. 优化初始化 =================
u = 0.46; v = 0.16; % 初始值
lambda_lm = 1e-1;   % 初始阻尼
loss_history = [];
iter_max = 100;

figure('Color', 'w', 'Position', [100 100 1000 400]);

for iter = 1:iter_max
    % --- A. 计算当前模型和残差 ---
    phase = 2*pi*d*((x_idx-cx)*u + (y_idx-cy)*v);
    O = exp(1i * phase);
    Z = O + R;
    r = abs(Z).^2 - I;
    
    % --- B. 显式定义损失函数 ---
    current_loss = 0.5 * sum(r(:).^2);
    loss_history(iter) = current_loss;
    
    % --- C. 计算雅可比和梯度 ---
    dO_du = 1i * 2*pi*d*(x_idx-cx) .* O;
    dO_dv = 1i * 2*pi*d*(y_idx-cy) .* O;
    J = [2*real(dO_du(:).*conj(Z(:))), 2*real(dO_dv(:).*conj(Z(:)))]; % 这里的 J 就是 Jacobian
    grad = J' * r(:); % 这就是梯度！指向 Loss 上升最快的方向
    
    % --- D. 尝试更新 (LM 步长) ---
    H = J' * J;
    step = -(H + lambda_lm * eye(2)) \ grad; % 沿着负梯度方向尝试
    
    u_trial = u + step(1);
    v_trial = v + step(2);
    
    % --- E. 验证更新是否有效 ---
    O_trial = exp(1i * 2*pi*d*((x_idx-cx)*u_trial + (y_idx-cy)*v_trial));
    loss_trial = 0.5 * sum( (abs(O_trial + R).^2 - I).^2, 'all' );
    
    if loss_trial < current_loss
        % 优化的方向是对的！
        u = u_trial; v = v_trial;
        lambda_lm = lambda_lm / 10; % 减小阻尼，加速前进
    else
        % 优化反了或跳过了！
        lambda_lm = lambda_lm * 10; % 增加阻尼，改走稳健的小步梯度下降
        continue; % 重新尝试本轮
    end
    
    % --- 可视化 ---
    if mod(iter, 2) == 0
        subplot(1,2,1); plot(loss_history, 'r-o', 'LineWidth', 1.5);
        xlabel('迭代次数'); ylabel('损失函数值 (Loss)'); title('收敛曲线'); grid on;
        subplot(1,2,2); imagesc(angle(O .* conj(O_true))); axis image; colorbar;
        title(['相位误差 (Iter ' num2str(iter) ')']);
        drawnow;
    end
    
    if norm(step) < 1e-7, break; end
end

fprintf('最终结果: u=%.4f (真值:%.4f), v=%.4f (真值:%.4f)\n', u, u_true, v, v_true);