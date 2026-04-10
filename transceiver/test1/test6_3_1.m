%% 完成了数据传输的验证功能
clc; clear; close all;

%% ================= System parameters =================
Ny = 5; Nz = 5;
K = Ny*Nz;
L = 2;
lambda = 1;
dy = lambda/2; dz = lambda/2;

b = ones(K,1);     % known reference wave

theta = [20, 40]*pi/180;
phi   = [-30, 10]*pi/180;
alpha = [1+0.5j; 0.8-0.3j];

V = zeros(K,L);
for l = 1:L
    V(:,l) = steering2D(theta(l),phi(l),Ny,Nz,dy,dz,lambda);
end
A = diag(alpha);

%% ================= QPSK data =================
T = 300;
S = qpsk_mod(L,T);

%% ================= Power measurements =================
I = zeros(K,T);
for t = 1:T
    I(:,t) = abs(V*A*S(:,t) + b).^2;
end

%% ================= GN symbol estimation =================
maxIter = 20;
err = zeros(L,1);

for t = 1:T

    % ---- initialization near zero (important)
    x = zeros(2*L,1);

    for iter = 1:maxIter
        [r,J] = residual_jacobian_symbol(x,V,A,b,I(:,t));
        dx = -(J.'*J)\(J.'*r);
        x = x + dx;

        if norm(dx) < 1e-6
            break;
        end
    end

    s_hat = x(1:L) + 1j*x(L+1:end);
    s_hat = qpsk_hard(s_hat);

    err = err + (s_hat ~= S(:,t));
end

SER = err/T;
disp('Symbol Error Rate per user:');
disp(SER);


%%
function S = qpsk_mod(L,T)
bits = randi([0 1],2*L,T);
S = (2*bits(1:2:end,:) - 1) ...
  + 1j*(2*bits(2:2:end,:) - 1);
end

function s = qpsk_hard(s)
s = sign(real(s)) + 1j*sign(imag(s));

% 防止 sign(0) = 0 的极端数值情况
s(real(s)==0) = 1 + 1j*sign(imag(s(real(s)==0)));
s(imag(s)==0) = sign(real(s(imag(s)==0))) + 1j;
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
    xp = x;
    xp(k) = xp(k) + eps;
    sp = xp(1:L) + 1j*xp(L+1:end);
    yp = V*A*sp + b;
    mup = abs(yp).^2;
    J(:,k) = (mup - mu)/eps;
end
end

function v = steering2D(theta,phi,Ny,Nz,dy,dz,lambda)

ky = 2*pi/lambda * sin(theta).*cos(phi);
kz = 2*pi/lambda * sin(theta).*sin(phi);

ay = exp(1j*ky*(0:Ny-1).'*dy);
az = exp(1j*kz*(0:Nz-1).'*dz);

v = kron(az,ay);
end
