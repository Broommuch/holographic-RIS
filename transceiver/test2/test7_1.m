%% 完美的doa估计和符号数据传输 在其基础上添加有无参考波的仿真对比图以及二维功率可视化的图
clc; clear; close all;

%% ================= RIS & system =================
Ny = 8; Nz = 8;
K = Ny*Nz;
L = 2;

lambda = 1;
dy = lambda/2; dz = lambda/2;

    %% 对已知参考波幅度，但未知参考波相位的测试，效果很差
    % phase = 2 * pi * rand(K, 1);    % 生成 [0, 2π) 范围内的随机相位
    % x = exp(1j * phase);            % 幅值为 1，相位为 phase 的复数向量
    % b = x;                  % reference wave
%%
b = ones(K,1);                  % reference wave
% b = zeros(K,1);                  % reference wave

SNRdB = 0:10:20;
MC = 10;

rng(1);
%% ================= True parameters =================
theta_true = [20, 60]*pi/180;
phi_true   = [-10, 10]*pi/180;
alpha_true = [1+0.5j; 0.8-0.3j];

%% ================= Pilot & data =================
Tp = 200;
Td = 200;

S_pilot = qpsk_mod(L,Tp);
S_data  = qpsk_mod(L,Td);

%% ================= Performance storage =================
MSE_theta = zeros(length(SNRdB),1);
MSE_phi   = zeros(length(SNRdB),1);
SER       = zeros(length(SNRdB),1);
crb_theta = zeros(length(SNRdB),1);
crb_phi   = zeros(length(SNRdB),1);

%% ================= SNR loop =================
for isnr = 1:length(SNRdB)

    snr = SNRdB(isnr);
    sigma2 = 10^(-snr/10);

    mse_theta_mc = 0;
    mse_phi_mc   = 0;
    ser_mc       = 0;
    total_sym    = 0;

    for mc = 1:MC

        %% ========== Pilot phase ==========
        I_p = zeros(K,Tp);
        for t = 1:Tp
            y = zeros(K,1);
            for l = 1:L
                v = steering2D(theta_true(l),phi_true(l),Ny,Nz,dy,dz,lambda);
                y = y + alpha_true(l)*v*S_pilot(l,t);
            end
            I_p(:,t) = abs(y + b).^2 + sqrt(sigma2)*randn(K,1);
            I_p_wb = abs(y).^2 + sqrt(sigma2)*randn(K,1);
        end

        % ---- GN initialization
        theta0 = theta_true + 5*pi/180*randn(1,L);
        phi0   = phi_true   + 5*pi/180*randn(1,L);
        alpha0 = alpha_true .* (1+0.2*(randn(L,1)+1j*randn(L,1)));

        x = pack_param(theta0,phi0,alpha0);

        % ---- GN iterations
        for iter = 1:25
            [r,J] = residual_jacobian_doa(x,I_p,S_pilot,abs(b),...
                                          Ny,Nz,dy,dz,lambda);
%             [r_wb,J_wb] = residual_jacobian_doa(x,I_p_wb,S_pilot,b,...
%                                           Ny,Nz,dy,dz,lambda);
%             dx = -(J.'*J)\(J.'*r);
            mu = 1e-3 * trace(J.'*J) / size(J,2);   % 自适应阻尼
            dx = -(J.'*J + mu*eye(size(J,2))) \ (J.'*r);

            x = x + dx;
            if norm(dx) < 1e-5, break; end
        end

        [theta_hat,phi_hat,alpha_hat] = unpack_param(x,L);

        mse_theta_mc = mse_theta_mc + mean((theta_hat-theta_true.').^2);
        mse_phi_mc   = mse_phi_mc   + mean((phi_hat-phi_true.').^2);

        %% ========== Construct estimated V, A ==========
        V_hat = zeros(K,L);
        for l = 1:L
            V_hat(:,l) = steering2D(theta_hat(l),phi_hat(l),...
                                    Ny,Nz,dy,dz,lambda);
        end
        A_hat = diag(alpha_hat);

        %% ========== Data phase ==========
        for t = 1:Td
            y = zeros(K,1);
            for l = 1:L
                v = steering2D(theta_true(l),phi_true(l),Ny,Nz,dy,dz,lambda);
                y = y + alpha_true(l)*v*S_data(l,t);
            end
            I_d = abs(y + b).^2 + sqrt(sigma2)*randn(K,1);

            s_hat = estimate_symbol_GN(V_hat,A_hat,abs(b),I_d);
            s_hat = qpsk_hard(s_hat);

            ser_mc = ser_mc + sum(s_hat ~= S_data(:,t));
            total_sym = total_sym + L;
        end
    end

    % === CRB computation for DOA ===
%     sigma2 = noise_var;   % 你已有
    F = (J' * J) / sigma2;
    CRB = diag(inv(F));
    
    crb_theta(isnr) = mean(CRB(1:L));
    crb_phi(isnr)   = mean(CRB(L+1:2*L));


    MSE_theta(isnr) = mse_theta_mc/MC;
    MSE_phi(isnr)   = mse_phi_mc/MC;
    SER(isnr)       = ser_mc/total_sym;

    fprintf('SNR = %d dB finished\n',snr);

    % 选一个中等 SNR 的 snapshot（例如最后一次 mc 的某一时刻）
    y_snapshot = I_d;   % 你在 data phase 里已有
    
    plot_ris_results( ...
        SNRdB, ...
        MSE_theta, MSE_phi, ...
        crb_theta, crb_phi, ...
        SER, ...
        Ny, Nz, ...
        y_snapshot, ...
        'GN_RIS');
    
    



end


% compare_refwave_doa( ...
%     SNRdB, MC, ...
%     theta_true, phi_true, alpha_true, ...
%     Ny, Nz, dy, dz, lambda, ...
%     S_pilot);

%% ================= Plot =================
% figure;
% semilogy(SNRdB,MSE_theta,'-o','LineWidth',1.5); hold on;
% semilogy(SNRdB,MSE_phi,'-s','LineWidth',1.5); hold on;
% semilogy(SNRdB, crb_theta, '--', 'LineWidth', 1.5); hold on;
% semilogy(SNRdB, crb_phi, '-*', 'LineWidth', 1.5);
% grid on;
% xlabel('SNR (dB)'); ylabel('MSE');
% legend('Elevation','Azimuth','Elevation_CRB','Azimuth_CRB');
% title('DOA MSE vs SNR');
% 
% figure;
% semilogy(SNRdB,SER,'-o','LineWidth',1.5);
% grid on;
% xlabel('SNR (dB)'); ylabel('SER');
% title('Symbol Error Rate vs SNR');





