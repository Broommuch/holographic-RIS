clc; clear; close all;

%% ================= 参数设置 =================
M = 16;
Nsym = 100;

Rs = 5e6;
sps = 20;
fs = Rs * sps;        % ✅ 正确采样率

fc = 3.5e9;           % 仅用于标注（不真实采样）

rolloff = 0.25;
span = 6;

bias = 1 + 0j;

rng(100);
%% ================= 阵列 & 用户参数 =================
K = 16;               % RIS单元数
L = 2;                % 用户数

d = 0.5;              % 阵元间距（单位：lambda）

theta = [20, -30] * pi/180;   % 方位角（弧度）
phi   = [10, 5] * pi/180;     % 俯仰角

alpha = [1, 0.8] .* exp(1j*[0, pi/4]);   % 复信道增益

%% ================= 构造导向矢量 =================
V = zeros(K, L);

for l = 1:L
    for k = 1:K
        % 简化ULA模型（只用方位角）
        V(k,l) = exp(1j * 2*pi * d * (k-1) * sin(theta(l)));
    end
end

H = V * diag(alpha);   % K x L 信道矩阵

%% ================= 生成每个用户信号 =================
rrc = rcosdesign(rolloff, span, sps, 'sqrt');

tx_bb_all = [];

for l = 1:L
    data = randi([0 M-1], Nsym, 1);
    symbols = qammod(data, M, 'UnitAveragePower', true);
    symbols_bias = symbols + bias;

    tx_bb = upfirdn(symbols_bias, rrc, sps, 1);

    tx_bb_all(:,l) = tx_bb;
end

Nt = size(tx_bb_all,1);
t = (0:Nt-1)/fs;

%% ================= 等效RF信号 =================
tx_rf_all = zeros(Nt, L);

for l = 1:L
    tx_rf_all(:,l) = real(tx_bb_all(:,l) .* exp(1j*2*pi*(fc/fs)*t.'));
end

%% ================= RIS接收信号 =================
g = zeros(Nt, K);

for k = 1:K
    for l = 1:L
        g(:,k) = g(:,k) + real( H(k,l) * tx_bb_all(:,l) .* exp(1j*2*pi*(fc/fs)*t.') );
    end
end

%% ================= 加参考信号 =================
A_ref = 1;
phi_ref = 0;

ref = A_ref * cos(2*pi*(fc/fs)*t + phi_ref);

for k = 1:K
    g(:,k) = g(:,k) + ref(:);
end

%% ================= 加噪声 =================
SNR_dB = 20;

signal_power = mean(g(:).^2);
noise_power = signal_power / (10^(SNR_dB/10));

noise = sqrt(noise_power) * randn(size(g));

g = g + noise;

%% ================= 平方律检测 =================
z = g.^2;

%% ================= 低通滤波 =================
% 截止频率设为基带带宽
z_lp = zeros(size(z));

for k = 1:K
    z_lp(:,k) = lowpass(z(:,k), Rs, fs);
end

%% ================= 能量输出可视化 =================
Z = mean(z_lp, 1).';   % K×1 能量向量




%% ================= 基带频谱 =================
Nfft = 4096;
f = linspace(-fs/2, fs/2, Nfft);

% 只看第一个用户的基带
S_bb = fftshift(abs(fft(tx_bb_all(:,1), Nfft)));



%% ================= RF频谱（平移到3.5GHz） =================
% 也只看第一个用户的射频频谱
S_rf = fftshift(abs(fft(tx_rf_all(:,1), Nfft)));

f_rf = f + fc;   % 平移频率轴

%% ================ 画图 ============
% ================ 能量输出可视化 =================
figure;
plot(Z, 'o-');
xlabel('RIS element index');
ylabel('Measured Energy');
title('Energy Measurement across RIS Array');
grid on;

% ================ 时域波形 =================
figure;
subplot(2,1,1);
plot(t(1:500)*1e6, tx_rf_all(1:500));
xlabel('Time (us)');
ylabel('Amplitude');
title('Equivalent RF Signal (3.5 GHz carrier)');
grid on;

subplot(2,1,2);
plot(t(1:500)*1e6, real(tx_bb(1:500)));
hold on;
plot(t(1:500)*1e6, imag(tx_bb(1:500)));
legend('I(t)','Q(t)');
title('Baseband Signal');
grid on;

% 基带频谱
figure;
plot(f/1e6, 20*log10(S_bb/max(S_bb)));
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
title('Baseband Spectrum');
grid on;

% 射频频谱
figure;
plot(f_rf/1e9, 20*log10(S_rf/max(S_rf)));
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
title('RF Spectrum centered at 3.5 GHz');
grid on;

figure;
% 多个单元波形
subplot(2,1,1);
plot(t(1:500)*1e6, g(1:500,1));
title('Received RF Signal at RIS Element 1');
xlabel('Time (us)');
grid on;
subplot(2,1,2);
plot(t(1:500)*1e6, g(1:500,2));
title('Received RF Signal at RIS Element 2');
xlabel('Time (us)');
grid on;