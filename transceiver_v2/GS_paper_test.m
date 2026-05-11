clc; clear; close all;

%% ================= 参数设置 =================
K = 50;        % 符号数
N = 150;       % 观测数，建议大于K
t0 = 300;      % GS迭代次数

rng(64);

%% ================= 生成QPSK符号 =================
bits = randi([0 1], 1, 2*K);
bits_reshape = reshape(bits, 2, []).';
symbols_idx = bi2de(bits_reshape, 'left-msb');

s_true = qammod(symbols_idx, 4, 'gray', ...
    'UnitAveragePower', true);

s_true = s_true(:);

%% ================= 构造已知矩阵A =================
A = randn(K, N) + 1i * randn(K, N);
A = A / sqrt(K);

%% ================= 构造已知参考信号 =================
% 不建议用QPSK对称方向，例如 pi/4
% 先用非对称相位参考，减少几何歧义
ref_amp = 1.2;
ref_phase = pi / 7;

b_ref = ref_amp * exp(1j * ref_phase) * ones(N, 1);

%% ================= 验证符号域GS恢复 =================
result = verify_symbol_domain_gs(s_true, b_ref, A, t0);