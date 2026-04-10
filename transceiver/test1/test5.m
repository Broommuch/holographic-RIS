% 全息接收机另外的版本，这次加上了载波和信道以及RIS的方位角的联合估计 

clc; clear; close all;

%% 基本参数
Ny = 8; Nz = 8;
dy = 0.5; dz = 0.5;
k = 2*pi;              % 归一化波数
K = Ny*Nz;

[y_idx, z_idx] = meshgrid(0:Ny-1, 0:Nz-1);
y = y_idx(:)*dy;
z = z_idx(:)*dz;

%% 真值
theta_true = 20/180*pi;
phi_true   = 40/180*pi;

ky = k*sin(theta_true)*cos(phi_true);
kz = k*sin(theta_true)*sin(phi_true);

alpha = 1.2*exp(1j*0.3);   % 已知
s_true = exp(1j*0.7);

b = 2*exp(1j*0.1);         % 已知参考波

%% 生成测量
noise_var = 1e-4;
I = zeros(K,1);
for i = 1:K
    h = alpha * s_true * exp(1j*(ky*y(i)+kz*z(i)));
    I(i) = abs(h + b)^2 + sqrt(noise_var)*randn;
end

%% GN 初值
theta0 = [0.8; 0.2; 15/180*pi; 30/180*pi]; 
% [Re{s}, Im{s}, theta, phi]

%% 高斯-牛顿
theta_hat = gauss_newton(I, y, z, b, alpha, k, theta0);

%% 结果
s_hat = theta_hat(1) + 1j*theta_hat(2);
theta_est = theta_hat(3);
phi_est   = theta_hat(4);

fprintf('True s = %.3f%+.3fj\n', real(s_true), imag(s_true));
fprintf('Est  s = %.3f%+.3fj\n', real(s_hat), imag(s_hat));
fprintf('Theta true/est = %.2f / %.2f deg\n', ...
        theta_true*180/pi, theta_est*180/pi);
fprintf('Phi   true/est = %.2f / %.2f deg\n', ...
        phi_true*180/pi, phi_est*180/pi);

function theta = gauss_newton(I, y, z, b, alpha, k, theta0)

theta = theta0;
maxIter = 50;

for iter = 1:maxIter
    [r, J] = residual_jacobian(theta, I, y, z, b, alpha, k);
    delta = -(J'*J)\(J'*r);
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
    
    % 对 Re{s}
    J(i,1) = dI_dh * real(alpha*phase);
    % 对 Im{s}
    J(i,2) = dI_dh * real(1j*alpha*phase);
    
    % 对 theta
    dky = k*cos(th)*cos(ph);
    dkz = k*cos(th)*sin(ph);
    dh_dth = alpha*s*1j*(dky*y(i)+dkz*z(i))*phase;
    J(i,3) = dI_dh * real(dh_dth);
    
    % 对 phi
    dky = -k*sin(th)*sin(ph);
    dkz =  k*sin(th)*cos(ph);
    dh_dph = alpha*s*1j*(dky*y(i)+dkz*z(i))*phase;
    J(i,4) = dI_dh * real(dh_dph);
end
end
