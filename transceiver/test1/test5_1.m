clc; clear; close all;

%% ================= 系统参数 =================
Ny = 8; Nz = 8;
dy = 0.5; dz = 0.5;
k = 2*pi;
K = Ny*Nz;

[y_idx, z_idx] = meshgrid(0:Ny-1, 0:Nz-1);
y = y_idx(:)*dy;
z = z_idx(:)*dz;

alpha = 1.2*exp(1j*0.3);      % 已知
b = 2*exp(1j*0.1);            % 已知参考波

theta_true = 20/180*pi;
phi_true   = 40/180*pi;

%% QPSK
QPSK = [1+1j, 1-1j, -1+1j, -1-1j]/sqrt(2);

%% 仿真参数
SNR_dB_list = -10:1:20;
Nsym = 500;

SER = zeros(length(SNR_dB_list),1);
theta_mse = zeros(length(SNR_dB_list),1);
phi_mse   = zeros(length(SNR_dB_list),1);

constellation_est = [];

%% ================= 主仿真 =================
for snr_idx = 1:length(SNR_dB_list)
    
    SNR_dB = SNR_dB_list(snr_idx);
    noise_var = 10^(-SNR_dB/10);
    
    sym_err = 0;
    theta_err = 0;
    phi_err = 0;
    
    for n = 1:Nsym
        
        %% 发送符号
        s_true = QPSK(randi(4));
        
        %% 生成测量
        I = zeros(K,1);
        ky = k*sin(theta_true)*cos(phi_true);
        kz = k*sin(theta_true)*sin(phi_true);
        
        for i = 1:K
            h = alpha*s_true*exp(1j*(ky*y(i)+kz*z(i)));
            I(i) = abs(h + b)^2 + sqrt(noise_var)*randn;
        end
        
        %% GN 初值
        theta0 = [0.8; 0.2; theta_true+0.05; phi_true+0.05];
        
        %% GN 估计
        theta_hat = gauss_newton(I, y, z, b, alpha, k, theta0);
        
        s_hat = theta_hat(1) + 1j*theta_hat(2);
        theta_est = theta_hat(3);
        phi_est   = theta_hat(4);
        
        %% QPSK 判决
        [~, idx] = min(abs(s_hat - QPSK.'));
        s_dec = QPSK(idx);
        
        %% 统计
        sym_err = sym_err + (s_dec ~= s_true);
        theta_err = theta_err + (theta_est-theta_true)^2;
        phi_err   = phi_err   + (phi_est-phi_true)^2;
        
        if SNR_dB == max(SNR_dB_list)
            constellation_est = [constellation_est; s_hat];
        end
    end
    
    SER(snr_idx) = sym_err / Nsym;
    theta_mse(snr_idx) = theta_err / Nsym;
    phi_mse(snr_idx)   = phi_err / Nsym;
    
    fprintf('SNR=%2d dB | SER=%.3e | theta RMSE=%.3e\n', ...
        SNR_dB, SER(snr_idx), sqrt(theta_mse(snr_idx)));
end


%%
figure;
semilogy(SNR_dB_list, SER, '-o','LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('SER');
title('Symbol Error Rate');

figure;
semilogy(SNR_dB_list, theta_mse*180^2/pi^2, '-s','LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Angle MSE (deg^2)');
title('Angle Estimation MSE');

figure;
plot(real(constellation_est), imag(constellation_est), 'r.');
hold on;
plot(real(QPSK), imag(QPSK), 'ko','MarkerSize',10,'LineWidth',2);
grid on; axis equal;
title('Estimated Constellation (High SNR)');

%%
function theta = gauss_newton(I, y, z, b, alpha, k, theta0)

theta = theta0;
maxIter = 30;

for iter = 1:maxIter
    [r, J] = residual_jacobian(theta, I, y, z, b, alpha, k);
    delta = -(J.'*J)\(J.'*r);
    theta = theta + delta;
    
    if norm(delta) < 1e-6
        break;
    end
end
end

function [r, J] = residual_jacobian(theta, I, y, z, b, alpha, k)

s = theta(1) + 1j*theta(2);
th = theta(3);
ph = theta(4);

ky = k*sin(th)*cos(ph);
kz = k*sin(th)*sin(ph);

K = length(I);
r = zeros(K,1);
J = zeros(K,4);

for i = 1:K
    phase = exp(1j*(ky*y(i)+kz*z(i)));
    h = alpha*s*phase;
    y_i = h + b;
    
    r(i) = I(i) - abs(y_i)^2;
    dI_dh = -2*real(conj(y_i));
    
    J(i,1) = dI_dh * real(alpha*phase);
    J(i,2) = dI_dh * real(1j*alpha*phase);
    
    dky = k*cos(th)*cos(ph);
    dkz = k*cos(th)*sin(ph);
    dh_dth = alpha*s*1j*(dky*y(i)+dkz*z(i))*phase;
    J(i,3) = dI_dh * real(dh_dth);
    
    dky = -k*sin(th)*sin(ph);
    dkz =  k*sin(th)*cos(ph);
    dh_dph = alpha*s*1j*(dky*y(i)+dkz*z(i))*phase;
    J(i,4) = dI_dh * real(dh_dph);
end
end
