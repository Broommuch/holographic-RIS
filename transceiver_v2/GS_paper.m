%% 这个脚本复现了里德堡原子接收机中的有偏GS算法，可用于后续的符号提取

clc; clear; close all;

%% 参数设置
K = 5;       % Dimension of symbol vector s
N = 36;      % Number of observations
t0 = 200;      % Number of iterations

% 生成随机测试数据
A = randn(K, N) + 1i*randn(K, N);  % Random complex channel matrix
s_true = exp(1i*2*pi*rand(K, 1));  % True symbol vector (complex unit modulus)
% b = randn(N, 1) + 1i*randn(N, 1);  % Random reference field
b = (0.5 + 1i*0.5)*ones(36,1);  % fixed reference field
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


