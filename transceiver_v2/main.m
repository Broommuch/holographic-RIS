%% 版本过时，看main_test1即可
clc; clear; close all;

%% ================= 参数设置 =================
M = 16;                  % 16-QAM
Nsym = 100;               % 符号数

Rs = 5e6;                % 符号率 5 MHz
sps = 20;                % 每符号采样点
fs = Rs * sps;           % 采样率 = 100 MHz

fc = 3.5e9;              % 载波频率（用于标注）

rolloff = 0.25;
span = 6;                % 滤波器在时域上截断的长度，单位是"符号周期”

bias = 1 + 0j;           % 非零均值（关键）

%% ================= QAM符号 =================
data = randi([0 M-1], Nsym, 1);
symbols = qammod(data, M, 'UnitAveragePower', true);

% 加偏置
symbols_bias = symbols + bias;

%% ================= 成形滤波 =================
rrc = rcosdesign(rolloff, span, sps, 'sqrt');

tx_bb = upfirdn(symbols_bias, rrc, sps, 1);

t = (0:length(tx_bb)-1)/fs;

%% ================= 等效RF信号 =================
% 注意：这里不直接生成3.5GHz cos，而是等效表示
tx_rf_equiv = real(tx_bb .* exp(1j*2*pi*(fc/fs)*t.'));

%% ================= 时域波形 =================
figure;
subplot(2,1,1);
plot(t(1:500)*1e6, tx_rf_equiv(1:500));
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

%% ================= 基带频谱 =================
Nfft = 4096;
f = linspace(-fs/2, fs/2, Nfft);

S_bb = fftshift(abs(fft(tx_bb, Nfft)));

figure;
plot(f/1e6, 20*log10(S_bb/max(S_bb)));
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
title('Baseband Spectrum');
grid on;

%% ================= RF频谱（平移到3.5GHz） =================
S_rf = fftshift(abs(fft(tx_rf_equiv, Nfft)));

f_rf = f + fc;   % 平移频率轴

figure;
plot(f_rf/1e9, 20*log10(S_rf/max(S_rf)));
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
title('RF Spectrum centered at 3.5 GHz');
grid on;

%% ================= 平方律检测 =================
z = tx_rf_equiv.^2;

% 低通滤波（提取能量项）
z_lp = lowpass(z, Rs, fs);

figure;
plot(t(1:500)*1e6, z_lp(1:500));
xlabel('Time (us)');
title('After Square-law + LPF');
grid on;