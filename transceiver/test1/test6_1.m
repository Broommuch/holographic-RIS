%% 这一版在doa估计的基础上进行了数据传输，但是误符号率很高，需继续修改
clc; clear; close all;

%% ================= RIS & system =================
Ny = 5; Nz = 5;
K = Ny*Nz;
L = 2;
lambda = 1;
dy = lambda/2; dz = lambda/2;

b = ones(K,1);                   % reference wave

%% ================= True parameters =================
theta_true = [20, 40]*pi/180;
phi_true   = [-30, 10]*pi/180;
alpha_true = [1+0.5j; 0.8-0.3j];

%% ================= Pilot =================
Tp = 200;
Sp = exp(1j*pi/2*randi([0 3],L,Tp));   % known pilots

Ip = zeros(K,Tp);
for t = 1:Tp
    y = zeros(K,1);
    for l = 1:L
        v = steering2D(theta_true(l),phi_true(l),Ny,Nz,dy,dz,lambda);
        y = y + alpha_true(l)*v*Sp(l,t);
    end
    Ip(:,t) = abs(y + b).^2;
end

%% ================= Stage I: DOA + A estimation =================
theta0 = theta_true + 5*pi/180*randn(1,L);
phi0   = phi_true   + 5*pi/180*randn(1,L);
alpha0 = alpha_true .* (1+0.2*(randn(L,1)+1j*randn(L,1)));

x = pack_param(theta0,phi0,alpha0);

maxIter = 30;
for iter = 1:maxIter
    [r,J] = residual_jacobian(x,Ip,Sp,b,Ny,Nz,dy,dz,lambda);
    dx = -(J'*J)\(J'*r);
    x = x + dx;
end

[theta_est,phi_est,alpha_est] = unpack_param(x,L);

%% ================= Build estimated V and A =================
V_hat = zeros(K,L);
for l = 1:L
    V_hat(:,l) = steering2D(theta_est(l),phi_est(l),Ny,Nz,dy,dz,lambda);
end
A_hat = diag(alpha_est);

disp('Estimated DOA (deg):')
disp([theta_est'*180/pi phi_est'*180/pi])
disp('Estimated alpha:')
disp(alpha_est)

%% ================= Stage II: Data transmission =================
Td = 200;
S_data = exp(1j*pi/2*randi([0 3],L,Td));

% SNR_dB = 20;
% snr = 10^(SNR_dB/10);
% sigma2 = norm(b)^2/(K*snr);


Id = zeros(K,Td);
for t = 1:Td
    y = zeros(K,1);
    for l = 1:L
        v = steering2D(theta_true(l),phi_true(l),Ny,Nz,dy,dz,lambda);
        y = y + alpha_true(l)*v*S_data(l,t);
    end
    Id(:,t) = abs(y + b).^2 ;
end

%% ================= Symbol estimation =================
S_hat = zeros(L,Td);
opts = optimoptions('fminunc','Display','off');

for t = 1:Td
%     fun = @(s) norm(abs(V_hat*A_hat*s + b).^2 - Id(:,t))^2;
%     s0 = exp(1j*pi/2*randi([0 3],L,1));
%     s_est = fminunc(fun,s0,opts);
% --- 将复符号拆成实变量 ---
    fun = @(x) norm( ...
        abs(V_hat*A_hat*(x(1:L) + 1j*x(L+1:end)) + b).^2 ...
        - Id(:,t) ).^2;
    
    % 实变量初值
    s0 = exp(1j*pi/2*randi([0 3],L,1));
    x0 = [real(s0); imag(s0)];
    
    % 优化
    x_est = fminunc(fun, x0, opts);
    
    % 还原复符号
    s_est = x_est(1:L) + 1j*x_est(L+1:end);
    
    % QPSK 判决
    S_hat(:,t) = qpsk_hard(s_est);
end

%% ================= Constellation =================
figure;
subplot(1,2,1)
plot(real(S_data(:)),imag(S_data(:)),'bo'); axis equal; grid on;
title('True QPSK symbols')

subplot(1,2,2)
plot(real(S_hat(:)),imag(S_hat(:)),'r*'); axis equal; grid on;
title('Estimated QPSK symbols')


%%
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

function s = qpsk_hard(x)
s = exp(1j*pi/2*round(angle(x)/(pi/2)));
end
