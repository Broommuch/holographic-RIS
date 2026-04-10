% 这个脚本想尝试完整估计信道参数和方向导向矢量最后进行数据传输
% 信道参数和方向导向矢量 估计准确，但是数据传输误码率很高，需继续修改
clc; clear; close all;

%% ================= System =================
Ny = 5; Nz = 5;
K = Ny*Nz;
L = 2;
lambda = 1;
dy = lambda/2; dz = lambda/2;
b = ones(K,1);

theta_true = [20, 40]*pi/180;
phi_true   = [-30, 10]*pi/180;
alpha_true = [1+0.5j; 0.8-0.3j];

Tp = 100;
Td = 100;
SNR_dB = 0:5:10;
MC = 10;

%% ================= Storage =================
mse_theta = zeros(length(SNR_dB),1);
mse_phi   = zeros(length(SNR_dB),1);
ser       = zeros(length(SNR_dB),1);
crb_theta = zeros(length(SNR_dB),1);
crb_phi   = zeros(length(SNR_dB),1);

%% ================= Monte Carlo =================
for is = 1:length(SNR_dB)
    snr = 10^(SNR_dB(is)/10);
    sigma2 = norm(b)^2/(K*snr);

    err_theta = 0; err_phi = 0;
    err_sym = 0; tot_sym = 0;

    FIM = zeros(4*L);

    for mc = 1:MC

        %% ----- Pilot -----
        Sp = exp(1j*pi/2*randi([0 3],L,Tp));
        Ip = zeros(K,Tp);

        for t = 1:Tp
            y = zeros(K,1);
            for l = 1:L
                v = steering2D(theta_true(l),phi_true(l),Ny,Nz,dy,dz,lambda);
                y = y + alpha_true(l)*v*Sp(l,t);
            end
            Ip(:,t) = abs(y+b).^2 + sqrt(sigma2)*randn(K,1);
        end

        %% ----- GN Estimation -----
        theta0 = theta_true + 5*pi/180*randn(1,L);
        phi0   = phi_true   + 5*pi/180*randn(1,L);
        alpha0 = alpha_true;

        x = pack_param(theta0,phi0,alpha0);

        for iter = 1:20
            [r,J] = residual_jacobian(x,Ip,Sp,b,Ny,Nz,dy,dz,lambda);
            x = x - (J'*J)\(J'*r);
        end

        [theta_hat,phi_hat,alpha_hat] = unpack_param(x,L);

        err_theta = err_theta + mean((theta_hat-theta_true').^2);
        err_phi   = err_phi   + mean((phi_hat-phi_true').^2);

        %% ----- CRB accumulation -----
        FIM = FIM + fim_conditional(theta_true,phi_true,alpha_true,...
            Sp,b,Ny,Nz,dy,dz,lambda,sigma2);

        %% ----- Data -----
        Sd = exp(1j*pi/2*randi([0 3],L,Td));
        Id = zeros(K,Td);

        Vhat = zeros(K,L);
        for l = 1:L
            Vhat(:,l) = steering2D(theta_hat(l),phi_hat(l),Ny,Nz,dy,dz,lambda);
        end
        Ahat = diag(alpha_hat);

        for t = 1:Td
            y = zeros(K,1);
            for l = 1:L
                v = steering2D(theta_true(l),phi_true(l),Ny,Nz,dy,dz,lambda);
                y = y + alpha_true(l)*v*Sd(l,t);
            end
            Id(:,t) = abs(y+b).^2 + sqrt(sigma2)*randn(K,1);

            % symbol estimation
            fun = @(x) norm(abs(Vhat*Ahat*(x(1:L)+1j*x(L+1:end))+b).^2 - Id(:,t))^2;
            x0 = [1;1;0;0];
            xest = fminunc(fun,x0,optimoptions('fminunc','Display','off'));
            sest = qpsk_hard(xest(1:L)+1j*xest(L+1:end));

            err_sym = err_sym + sum(sest~=Sd(:,t));
            tot_sym = tot_sym + L;
        end
    end

    mse_theta(is) = err_theta/MC;
    mse_phi(is)   = err_phi/MC;
    ser(is)       = err_sym/tot_sym;

    CRB = inv(FIM/MC);
    crb_theta(is) = mean(diag(CRB(1:L,1:L)));
    crb_phi(is)   = mean(diag(CRB(L+1:2*L,L+1:2*L)));

    fprintf('SNR=%2d dB \n', ...
    snr);
end

%% ================= Plots =================
figure;
semilogy(SNR_dB,mse_theta,'o-',SNR_dB,crb_theta,'--','LineWidth',1.5);
grid on; xlabel('SNR (dB)'); ylabel('MSE');
legend('MSE \theta','CRB \theta');

figure;
semilogy(SNR_dB,mse_phi,'o-',SNR_dB,crb_phi,'--','LineWidth',1.5);
grid on; xlabel('SNR (dB)'); ylabel('MSE');
legend('MSE \phi','CRB \phi');

figure;
semilogy(SNR_dB,ser,'o-','LineWidth',1.5);
grid on; xlabel('SNR (dB)'); ylabel('SER');


%%
function F = fim_conditional(theta,phi,alpha,S,b,Ny,Nz,dy,dz,lambda,sigma2)

L = length(theta); T = size(S,2);
P = 4*L;
F = zeros(P);

eps = 1e-6;
x = pack_param(theta,phi,alpha);

for t = 1:T
    mu0 = intensity(x,S(:,t),b,Ny,Nz,dy,dz,lambda);
    G = zeros(length(mu0),P);

    for p = 1:P
        xp = x; xp(p) = xp(p)+eps;
        mup = intensity(xp,S(:,t),b,Ny,Nz,dy,dz,lambda);
        G(:,p) = (mup-mu0)/eps;
    end

    F = F + (G'*G)/sigma2;
end
end

function mu = intensity(x,s,b,Ny,Nz,dy,dz,lambda)
[theta,phi,alpha] = unpack_param(x,length(s));
K = Ny*Nz;
y = zeros(K,1);
for l = 1:length(s)
    v = steering2D(theta(l),phi(l),Ny,Nz,dy,dz,lambda);
    y = y + alpha(l)*v*s(l);
end
mu = abs(y+b).^2;
end


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

function [r,J] = residual_jacobian(x,I,S,b,Ny,Nz,dy,dz,lambda)

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

function [theta,phi,alpha] = unpack_param(x,L)
theta = x(1:L);
phi   = x(L+1:2*L);
alpha = x(2*L+1:3*L) + 1j*x(3*L+1:4*L);
end

function s = qpsk_hard(x)
s = exp(1j*pi/2*round(angle(x)/(pi/2)));
end