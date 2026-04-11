clc; clear; close all;

%% ================= 参数设置 =================
M = 16;                    % 16-QAM
Nsym = 100;                 % 符号数
Rs = 1e6;                  % 符号率 1 MHz
sps = 20;                  % 每符号采样点（过采样）
fs = Rs * sps;             % 采样率（必须满足 fs >> fc 的条件后面会统一）

fc = 10e6;                 % 载波频率 10 MHz（不要太高，否则仿真压力大）

rolloff = 0.25;            % RRC滚降系数
span = 6;                  % 滤波器长度（单位：符号）

bias = 1 + 0j;             % 偏置（关键！！非零均值）

%% ================= 生成QAM符号 =================
data = randi([0 M-1], Nsym, 1);
symbols = qammod(data, M, 'UnitAveragePower', true);

%% ================= 加偏置 =================
symbols_bias = symbols + bias;

%% ================= 成形滤波 =================
rrc = rcosdesign(rolloff, span, sps, 'sqrt');

% 上采样 + 滤波
tx_baseband = upfirdn(symbols_bias, rrc, sps, 1);

t = (0:length(tx_baseband)-1)/fs;

%% ================= 上变频到RF =================
% 复包络
s_complex = tx_baseband;

% RF信号（实信号）
tx_rf = real(s_complex .* exp(1j*2*pi*fc*t.'));

%% ================= 时域波形 =================
figure;
subplot(2,1,1);
plot(t(1:500)*1e6, tx_rf(1:500));
xlabel('Time (us)');
ylabel('Amplitude');
title('RF Signal (Time Domain)');
grid on;

subplot(2,1,2);
plot(t(1:500)*1e6, real(s_complex(1:500)));
hold on;
plot(t(1:500)*1e6, imag(s_complex(1:500)));
legend('I(t)','Q(t)');
title('Baseband I/Q');
grid on;

%% ================= 频谱分析 =================
Nfft = 4096;
f = linspace(-fs/2, fs/2, Nfft);

S = fftshift(abs(fft(tx_rf, Nfft)));

figure;
plot(f/1e6, 20*log10(S/max(S)));
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
title('RF Spectrum');
grid on;

%% ================= 基带频谱 =================
S_bb = fftshift(abs(fft(s_complex, Nfft)));

figure;
plot(f/1e6, 20*log10(S_bb/max(S_bb)));
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
title('Baseband Spectrum');
grid on;