%% 这个脚本简单验证了数据传输问题，但是qpsk习惯不对，需调整
clc; clear; close all;

%% ============ Known system ============
Ny = 5; Nz = 5;
K = Ny*Nz;
L = 2;
lambda = 1;
dy = lambda/2; dz = lambda/2;

b = ones(K,1);

theta = [20, 40]*pi/180;
phi   = [-30, 10]*pi/180;
alpha = [1+0.5j; 0.8-0.3j];

V = zeros(K,L);
for l = 1:L
    V(:,l) = steering2D(theta(l),phi(l),Ny,Nz,dy,dz,lambda);
end
A = diag(alpha);

%% ============ Data ============
T = 200;
S = exp(1j*pi/2*randi([0 3],L,T));

I = zeros(K,T);
for t = 1:T
    I(:,t) = abs(V*A*S(:,t) + b).^2;
end

%% ============ GN symbol estimation ============
maxIter = 20;
SER = zeros(L,1);
err = zeros(L,1);

for t = 1:T

    % ---- initialization (important!)
    x = zeros(2*L,1);  % near origin is OK due to b

    for iter = 1:maxIter
        [r,J] = residual_jacobian_symbol(x,V,A,b,I(:,t));
        dx = -(J'*J)\(J'*r);
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



function s = qpsk_hard(s)
tol = 1e-2;
if real(s)<tol
    real(s) = 0;
end
if imag(s)<tol
    imag(s) = 0;
end
s = sign(real(s)) + 1j*sign(imag(s));
s = s/sqrt(2);
end


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