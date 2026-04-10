%% 这一版为了验证可行性只做了doa估计
clc; clear; close all;

%% ================= RIS & system =================
Ny = 5; Nz = 5;
K = Ny*Nz;
L = 2;
lambda = 1;
dy = lambda/2; dz = lambda/2;

T = 200;                         % pilot length
b = ones(K,1);                   % reference wave

%% ================= True parameters =================
theta_true = [20, 60]*pi/180;
phi_true   = [-10, 10]*pi/180;
alpha_true = [1+0.5j; 0.8-0.3j];

%% ================= Pilot =================
S = exp(1j*pi/2*randi([0 3],L,T));   % known QPSK pilots

%% ================= Generate measurements =================
I = zeros(K,T);
for t = 1:T
    y = zeros(K,1);
    for l = 1:L
        v = steering2D(theta_true(l),phi_true(l),Ny,Nz,dy,dz,lambda);
        y = y + alpha_true(l)*v*S(l,t);
    end
    I(:,t) = abs(y + b).^2;
end

%% ================= Initial guess =================
theta0 = (theta_true + 5*pi/180.*randn(1,L));
phi0   = (phi_true   + 5*pi/180.*randn(1,L));
alpha0 = alpha_true .* (1+0.2*(randn(L,1)+1j*randn(L,1)));

x0 = pack_param(theta0,phi0,alpha0);

%% ================= Gauss-Newton =================
maxIter = 30;
for iter = 1:maxIter
    [r,J] = residual_jacobian(x0,I,S,b,Ny,Nz,dy,dz,lambda);
    dx = -(J'*J)\(J'*r);
    x0 = x0 + dx;
    fprintf('Iter %d, residual = %.3e\n',iter,norm(r));
end

%% ================= Results =================
[theta_est,phi_est,alpha_est] = unpack_param(x0,L);

disp('True DOA (deg):')
disp([theta_true*180/pi phi_true*180/pi])
disp('Estimated DOA (deg):')
disp([theta_est'*180/pi phi_est'*180/pi])

disp('True alpha:')
disp(alpha_true)
disp('Estimated alpha:')
disp(alpha_est)


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
    V = zeros(K,L);
    for l = 1:L
        V(:,l) = steering2D(theta(l),phi(l),Ny,Nz,dy,dz,lambda);
        y = y + alpha(l)*V(:,l)*S(l,t);
    end
    e = abs(y+b).^2 - I(:,t);
    r(idx:idx+K-1) = e;

    % Jacobian (numerical)
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
