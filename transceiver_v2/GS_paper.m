clc; clear; close all;

%% 参数设置
K = 5;       % Dimension of symbol vector s
N = 36;      % Number of observations
t0 = 50;      % Number of iterations

% 生成随机测试数据
A = randn(K, N) + 1i*randn(K, N);  % Random complex channel matrix
s_true = exp(1i*2*pi*rand(K, 1));  % True symbol vector (complex unit modulus)
b = randn(N, 1) + 1i*randn(N, 1);  % Random reference field
z = abs(A' * s_true + b);          % Observed magnitudes

% 运行算法
s_est = biased_gs_algorithm(z, A, b, t0);

% 显示结果
disp('True symbols:');
disp(s_true.');
disp('Estimated symbols:');
disp(s_est.');
error = norm(s_est - s_true)/norm(s_true);
disp(['Relative error: ', num2str(error)]);


