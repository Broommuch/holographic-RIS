%% 这个脚本尝试直接用ghz级别的采样率直采射频信号，但是没有写完之后再看看

clc; clear; close all;

%% ================= 参数 =================
M = 16;
Nsym = 100;              % 减少符号数

Rs = 10e6;             % 提高符号率（减少采样点）
sps = 10;              
fs = Rs * sps;         % 100 MHz（仍远小于10 GHz）

fc = 3.5e9;            % 真正RF

bias = 1;

rolloff = 0.25;
span = 4;

%% ================= 信号 =================
data = randi([0 M-1], Nsym, 1);
symbols = qammod(data, M, 'UnitAveragePower', true);
symbols_bias = symbols + bias;

rrc = rcosdesign(rolloff, span, sps, 'sqrt');
tx_baseband = upfirdn(symbols_bias, rrc, sps, 1);

t = (0:length(tx_baseband)-1)/fs;

%% ================= RF =================
tx_rf = real(tx_baseband .* exp(1j*2*pi*fc*t.'));

%% ================= 时域 =================
figure;
plot(t(1:200)*1e6, tx_rf(1:200));
title('3.5 GHz RF (Aliased View)');
grid on;

%% ================= 频谱 =================
Nfft = 2048;
S = fftshift(abs(fft(tx_rf, Nfft)));
f = linspace(-fs/2, fs/2, Nfft);

figure;
plot(f/1e6, 20*log10(S/max(S)));
title('Spectrum (Aliased)');
grid on;