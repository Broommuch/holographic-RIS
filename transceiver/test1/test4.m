%% 结果不可行

clc; clear; close all;

%% 参数设置
N = 64;                 % RIS 单元数
SNR_dB = 25;
K = 3;                  % 相移数

theta = [0, 2*pi/3, 4*pi/3];
b_amp = 2;

%% QPSK 符号
QPSK = [1+1j, 1-1j, -1+1j, -1-1j]/sqrt(2);
s_true = QPSK(randi(4));

%% 信道
h_true = (randn(N,1)+1j*randn(N,1))/sqrt(2);

%% 测量
I = zeros(N,K);
noise_var = 10^(-SNR_dB/10);

for k = 1:K
    b = b_amp * exp(1j*theta(k));
    for i = 1:N
        y = h_true(i)*s_true + b;
        I(i,k) = abs(y)^2 + sqrt(noise_var)*randn;
    end
end

%% Step 1: 恢复 z_i = h_i s
z_est = recover_z_holographic(I, b_amp, theta);

%% Step 2: QPSK 判决
[s_hat, h_hat] = joint_qpsk_channel_est(z_est, QPSK);

%% 性能评估
fprintf('True s  = %.2f%+.2fj\n', real(s_true), imag(s_true));
fprintf('Est  s  = %.2f%+.2fj\n', real(s_hat), imag(s_hat));

figure;
plot(real(h_true), imag(h_true), 'bo'); hold on;
plot(real(h_hat), imag(h_hat), 'r+');
legend('True h_i', 'Estimated h_i');
axis equal; grid on;
title('Channel Estimation via Holographic RIS');

function z = recover_z_holographic(I, b_amp, theta)

N = size(I,1);

A = [ ...
    1, cos(theta(1)), sin(theta(1));
    1, cos(theta(2)), sin(theta(2));
    1, cos(theta(3)), sin(theta(3)) ];

z = zeros(N,1);

for i = 1:N
    y = I(i,:).' - b_amp^2;
    x = A \ y;
    z(i) = (x(2) + 1j*x(3)) / (2*b_amp);
end
end

function [s_hat, h_hat] = joint_qpsk_channel_est(z, QPSK)

M = length(QPSK);
metric = zeros(M,1);

for m = 1:M
    s = QPSK(m);
    h = z / s;
    metric(m) = sum(abs(h).^2);
end

[~, idx] = min(metric);
s_hat = QPSK(idx);
h_hat = z / s_hat;

end
