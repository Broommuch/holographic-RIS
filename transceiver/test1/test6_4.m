%% 报了内存错误的问题，还需要再改
clc; clear; close all;

%% =================== System parameters ===================
Ny = 5; Nz = 5;
K = Ny*Nz;
L = 2;

lambda = 1;
dy = lambda/2; dz = lambda/2;

b = ones(K,1);            % known reference wave
SNRdB = 0:5:30;
MC = 10;

%% True channel parameters
theta_true = [20, 40]*pi/180;
phi_true   = [-30, 10]*pi/180;
alpha_true = [1+0.5j; 0.8-0.3j];

V_true = zeros(K,L);
for l = 1:L
    V_true(:,l) = steering2D(theta_true(l),phi_true(l),Ny,Nz,dy,dz,lambda);
end
A_true = diag(alpha_true);

%% Pilot & data parameters
Tp = 20;          % pilot length
Td = 200;         % data length

S_pilot = qpsk_mod(L,Tp);
S_data  = qpsk_mod(L,Td);

%% Performance storage
MSE_theta = zeros(length(SNRdB),1);
MSE_phi   = zeros(length(SNRdB),1);
SER       = zeros(length(SNRdB),1);

%% =================== SNR loop ===================
for isnr = 1:length(SNRdB)

    snr = SNRdB(isnr);
    sigma2 = 10^(-snr/10);

    mse_theta_mc = 0;
    mse_phi_mc   = 0;
    ser_mc       = 0;
    total_sym    = 0;

    for mc = 1:MC

        %% ============ Pilot phase ============
        I_pilot = zeros(K,Tp);
        for t = 1:Tp
            y = V_true*A_true*S_pilot(:,t) + b;
            I_pilot(:,t) = abs(y).^2 + sqrt(sigma2)*randn(K,1);
        end

        % --- GN DOA + alpha estimation
        x0 = zeros(4*L,1);  % [theta phi Re(alpha) Im(alpha)]
        x0(1:2*L) = [theta_true(:); phi_true(:)] ...
                    + 5*pi/180*randn(2*L,1);

        opts.maxIter = 20;
        x = x0;

        for iter = 1:opts.maxIter
            [r,J] = residual_jacobian_doa(x,S_pilot,I_pilot,...
                                          Ny,Nz,dy,dz,lambda,b);
            dx = -(J.'*J)\(J.'*r);
            x = x + dx;
            if norm(dx) < 1e-5, break; end
        end

        theta_hat = x(1:L);
        phi_hat   = x(L+1:2*L);
        alpha_hat = x(2*L+1:3*L) + 1j*x(3*L+1:end);

        % --- construct estimated V and A
        V_hat = zeros(K,L);
        for l = 1:L
            V_hat(:,l) = steering2D(theta_hat(l),phi_hat(l),...
                                    Ny,Nz,dy,dz,lambda);
        end
        A_hat = diag(alpha_hat);

        mse_theta_mc = mse_theta_mc + mean((theta_hat-theta_true.').^2);
        mse_phi_mc   = mse_phi_mc   + mean((phi_hat-phi_true.').^2);

        %% ============ Data phase ============
        for t = 1:Td
            y = V_true*A_true*S_data(:,t) + b;
            I = abs(y).^2 + sqrt(sigma2)*randn(K,1);

            s_hat = estimate_symbol_GN(V_hat,A_hat,b,I);
            s_hat = qpsk_hard(s_hat);

            ser_mc = ser_mc + sum(s_hat ~= S_data(:,t));
            total_sym = total_sym + L;
        end
    end

    MSE_theta(isnr) = mse_theta_mc/(MC);
    MSE_phi(isnr)   = mse_phi_mc/(MC);
    SER(isnr)       = ser_mc/total_sym;

    fprintf('SNR = %d dB done\n',snr);
end

%% =================== Plots ===================
figure;
semilogy(SNRdB,MSE_theta,'-o','LineWidth',1.5); hold on;
semilogy(SNRdB,MSE_phi,'-s','LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('MSE');
legend('Elevation MSE','Azimuth MSE');
title('DOA estimation performance');

figure;
semilogy(SNRdB,SER,'-o','LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('SER');
title('QPSK symbol error rate');

%%
function s_hat = estimate_symbol_GN(V,A,b,I)

L = size(V,2);
x = zeros(2*L,1);
maxIter = 20;

for iter = 1:maxIter
    [r,J] = residual_jacobian_symbol(x,V,A,b,I);
    dx = -(J.'*J)\(J.'*r);
    x = x + dx;
    if norm(dx) < 1e-6, break; end
end

s_hat = x(1:L) + 1j*x(L+1:end);
end

function [r,J] = residual_jacobian_doa(x,S,I,Ny,Nz,dy,dz,lambda,b)

[L,Tp] = size(S);
theta = x(1:L);
phi   = x(L+1:2*L);
alpha = x(2*L+1:3*L) + 1j*x(3*L+1:end);

K = Ny*Nz;
V = zeros(K,L);
for l = 1:L
    V(:,l) = steering2D(theta(l),phi(l),Ny,Nz,dy,dz,lambda);
end
A = diag(alpha);

r = [];
for t = 1:Tp
    y = V*A*S(:,t) + b;
    r = [r; abs(y).^2 - I(:,t)];
end

eps = 1e-6;
J = zeros(length(r),length(x));

for k = 1:length(x)
    xp = x; xp(k) = xp(k) + eps;
    rp = residual_jacobian_doa(xp,S,I,Ny,Nz,dy,dz,lambda,b);
    J(:,k) = (rp - r)/eps;
end
end

function S = qpsk_mod(L,T)
bits = randi([0 1],2*L,T);
S = (2*bits(1:2:end,:) - 1) ...
  + 1j*(2*bits(2:2:end,:) - 1);
end

function s = qpsk_hard(s)
s = sign(real(s)) + 1j*sign(imag(s));
s(real(s)==0) = 1 + 1j*sign(imag(s(real(s)==0)));
s(imag(s)==0) = sign(real(s(imag(s)==0))) + 1j;
end

function v = steering2D(theta,phi,Ny,Nz,dy,dz,lambda)
ky = 2*pi/lambda*sin(theta).*cos(phi);
kz = 2*pi/lambda*sin(theta).*sin(phi);
ay = exp(1j*ky*(0:Ny-1).'*dy);
az = exp(1j*kz*(0:Nz-1).'*dz);
v = kron(az,ay);
end
