%% 完美的doa估计和符号数据传输
clc; clear; close all;

%% ================= RIS & system =================
Ny = 5; Nz = 5;
K = Ny*Nz;
L = 2;

lambda = 1;
dy = lambda/2; dz = lambda/2;

b = ones(K,1);                  % reference wave

SNRdB = 0:10:30;
MC = 10;

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
        end

        % ---- GN initialization
        theta0 = theta_true + 5*pi/180*randn(1,L);
        phi0   = phi_true   + 5*pi/180*randn(1,L);
        alpha0 = alpha_true .* (1+0.2*(randn(L,1)+1j*randn(L,1)));

        x = pack_param(theta0,phi0,alpha0);

        % ---- GN iterations
        for iter = 1:25
            [r,J] = residual_jacobian_doa(x,I_p,S_pilot,b,...
                                          Ny,Nz,dy,dz,lambda);
            dx = -(J.'*J)\(J.'*r);
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

            s_hat = estimate_symbol_GN(V_hat,A_hat,b,I_d);
            s_hat = qpsk_hard(s_hat);

            ser_mc = ser_mc + sum(s_hat ~= S_data(:,t));
            total_sym = total_sym + L;
        end
    end

    MSE_theta(isnr) = mse_theta_mc/MC;
    MSE_phi(isnr)   = mse_phi_mc/MC;
    SER(isnr)       = ser_mc/total_sym;

    fprintf('SNR = %d dB finished\n',snr);
end

%% ================= Plot =================
figure;
semilogy(SNRdB,MSE_theta,'-o','LineWidth',1.5); hold on;
semilogy(SNRdB,MSE_phi,'-s','LineWidth',1.5);
grid on;
xlabel('SNR (dB)'); ylabel('MSE');
legend('Elevation','Azimuth');
title('DOA MSE vs SNR');

figure;
semilogy(SNRdB,SER,'-o','LineWidth',1.5);
grid on;
xlabel('SNR (dB)'); ylabel('SER');
title('Symbol Error Rate vs SNR');


function v = steering2D(theta,phi,Ny,Nz,dy,dz,lambda)
ky = 2*pi/lambda*sin(theta)*cos(phi);
kz = 2*pi/lambda*sin(theta)*sin(phi);
ay = exp(1j*ky*(0:Ny-1).'*dy);
az = exp(1j*kz*(0:Nz-1).'*dz);
v = kron(az,ay);
v = v / norm(v);
end

function x = pack_param(theta,phi,alpha)
x = [theta(:); phi(:); real(alpha); imag(alpha)];
end

function [theta,phi,alpha] = unpack_param(x,L)
theta = x(1:L);
phi   = x(L+1:2*L);
alpha = x(2*L+1:3*L) + 1j*x(3*L+1:4*L);
end

function [r,J] = residual_jacobian_doa(x,I,S,b,Ny,Nz,dy,dz,lambda)

L = size(S,1); T = size(S,2);
[theta,phi,alpha] = unpack_param(x,L);
K = Ny*Nz;

r = zeros(K*T,1);
J = zeros(K*T,length(x));
idx = 1;

for t = 1:T
    y = zeros(K,1);
    for l = 1:L
        v = steering2D(theta(l),phi(l),Ny,Nz,dy,dz,lambda);
        y = y + alpha(l)*v*S(l,t);
    end
    e = abs(y+b).^2 - I(:,t);
    r(idx:idx+K-1) = e;

    eps = 1e-6;
    for p = 1:length(x)
        xp = x; xp(p) = xp(p) + eps;
        yp = zeros(K,1);
        [thp,php,alp] = unpack_param(xp,L);
        for l = 1:L
            vp = steering2D(thp(l),php(l),Ny,Nz,dy,dz,lambda);
            yp = yp + alp(l)*vp*S(l,t);
        end
        ep = abs(yp+b).^2 - I(:,t);
        J(idx:idx+K-1,p) = (ep - e)/eps;
    end
    idx = idx + K;
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

function s_hat = estimate_symbol_GN(V,A,b,I)
L = size(V,2);
x = zeros(2*L,1);
for iter = 1:20
    [r,J] = residual_jacobian_symbol(x,V,A,b,I);
    dx = -(J.'*J)\(J.'*r);
    x = x + dx;
    if norm(dx) < 1e-6, break; end
end
s_hat = x(1:L) + 1j*x(L+1:end);
end

function [r,J] = residual_jacobian_symbol(x,V,A,b,I)
L = length(x)/2;
s = x(1:L) + 1j*x(L+1:end);
y = V*A*s + b;
mu = abs(y).^2;
r = mu - I;

eps = 1e-6;
K = length(I);
J = zeros(K,2*L);

for k = 1:2*L
    xp = x; xp(k) = xp(k) + eps;
    sp = xp(1:L) + 1j*xp(L+1:end);
    yp = V*A*sp + b;
    mup = abs(yp).^2;
    J(:,k) = (mup - mu)/eps;
end
end
