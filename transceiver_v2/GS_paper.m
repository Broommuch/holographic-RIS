clc; clear; close all;

%% 参数设置
K = 10;       % Dimension of symbol vector s
N = 100;      % Number of observations
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


function s_est = biased_gs_algorithm(z, A, b, t0)
% Biased GS Algorithm for Phase Retrieval
% Inputs:
%   z: Equivalent received signal (N x 1 vector)
%   A: Equivalent channel matrix (K x N matrix)
%   b: Reference field (N x 1 vector)
%   t0: Total number of iterations (scalar)
% Output:
%   s_est: Estimated symbol vector (K x 1)

%% Step 1: Construct augmented observation matrix
A_H = A';               % A^H: conjugate transpose of A (N x K)
aug_mat = [A_H, b];     % Augmented matrix [A^H, b] (N x (K+1))
bar_A = aug_mat';       % \bar{A}: (K+1) x N matrix

%% Step 2: Find principal eigenvector of M
N = length(z);          % Number of observations
M = zeros(size(bar_A, 1)); % Initialize M ((K+1) x (K+1))

for n = 1:N
    a_n = bar_A(:, n);              % nth column of \bar{A}
    M = M + z(n) * (a_n * a_n');    % M = Σ z_n * (a_n * a_n^H)
end

[V, D] = eig(M);                   % Eigen decomposition
[~, idx] = max(diag(D));           % Index of largest eigenvalue
v = V(:, idx);                     % Principal eigenvector (K+1 x 1)

%% Step 3: Set initial vectors
temp = v' * bar_A;                 % 1 x N vector
numerator = abs(temp) * z;         % Scalar: |v^H \bar{A}| z
denominator = norm(bar_A' * v)^2;  % Scalar: ||\bar{A}^H v||_2^2
r_bar = numerator / denominator;   % Scalar
s_bar0 = r_bar * v;                % Initial estimate (K+1 x 1)

%% Step 4: Initialize s0
phase_shift = exp(-1i * angle(s_bar0(end)));  % Phase correction factor
s0 = phase_shift * s_bar0(1:end-1);           % Initial symbol estimate (K x 1)

%% Iteration loop
s_prev = s0;  % Initialize previous estimate

for t = 1:t0
    % Step 6: Update theta^t
    y = A' * s_prev + b;          % A^H s^{t-1} + b (N x 1)
    theta_t = angle(y);           % Element-wise phase angles
    
    % Step 7: Update s^t
    modulated = z .* exp(1i * theta_t);  % z ∘ e^{iθ^t}
    residual = modulated - b;            % z∘e^{iθ^t} - b
    s_t = (A * A') \ (A * residual);     % (AA^H)^{-1}A(residual)
    
    s_prev = s_t;  % Update for next iteration
end

s_est = s_prev;  % Final estimate after t0 iterations
end