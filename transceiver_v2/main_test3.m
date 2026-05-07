%% 这个脚本也几乎不能用，得重新调整

%% ================= 等效RF频率设置 =================
% 注意：fs = 100 MHz，不能真实采样 3.5 GHz 射频信号
% 因此这里用一个等效中频 fIF 来表示同频叠加过程
fIF = 10e6;     % 等效RF/IF频率，必须小于 fs/2
% fc = 3.5e9 仅作为物理工作频率标注，不直接进入离散采样

%% ================= 构造RIS接收到的未知RF信号 =================
g_sig = zeros(Nt, K);   % 不含参考信号的RIS接收RF信号

for k = 1:K
    for l = 1:L
        % 第 l 个用户到第 k 个RIS单元的复基带等效信号
        rx_bb_kl = H(k,l) * tx_bb_all(:,l);

        % 上变到等效RF/IF
        g_sig(:,k) = g_sig(:,k) + real(rx_bb_kl .* exp(1j*2*pi*fIF*t.'));
    end
end

%% ================= 构造已知同频参考信号 =================
A_ref = 1;          % 参考信号幅度
phi_ref = 0;        % 参考信号初始相位

% 所有RIS单元共用同一个参考信号
ref = A_ref * cos(2*pi*fIF*t.' + phi_ref);   % Nt × 1

% 如果希望每个RIS单元参考相位不同，可以用下面这种形式
% phi_ref_k = zeros(1,K);
% ref_all = zeros(Nt,K);
% for k = 1:K
%     ref_all(:,k) = A_ref * cos(2*pi*fIF*t.' + phi_ref_k(k));
% end

%% ================= 未知信号 + 已知参考信号相干叠加 =================
g = zeros(Nt, K);

for k = 1:K
    g(:,k) = g_sig(:,k) + ref;
end

%% ================= 加噪声 =================
SNR_dB = 20;

signal_power = mean(g(:).^2);
noise_power = signal_power / (10^(SNR_dB/10));

noise = sqrt(noise_power) * randn(size(g));
g = g + noise;

%% ================= 平方律检测 =================
% 平方律检波：输出中包含
% g^2 = signal^2 + ref^2 + 2 signal * ref
% 其中 2 signal * ref 是你想利用的相干交叉项
z = g.^2;

%% ================= 低通滤波，提取能量包络 =================
z_lp = zeros(size(z));

for k = 1:K
    % 低通截止频率可以取符号率附近
    % 若希望更平滑，可以改成 Rs/2 或 0.8*Rs
    z_lp(:,k) = lowpass(z(:,k), Rs, fs);
end

%% ================= 去除参考信号自身的直流能量项，可选 =================
% 参考信号平方后低通会产生 A_ref^2 / 2 的直流项
% 如果只关心未知信号导致的能量变化，可以减去该常数项
z_lp_remove_ref = z_lp - A_ref^2/2;

%% ================= 按符号周期积分，得到符号化能量检测结果 =================
% RRC滤波器群时延
gd = span * sps / 2;

% 每个符号中心对应的位置
sym_center = gd + 1 + (0:Nsym-1) * sps;

% 每个符号能量积分窗口
win_len = sps;
half_win = floor(win_len/2);

E_sym = zeros(Nsym, K);              % 未扣除参考直流项的符号能量
E_sym_remove_ref = zeros(Nsym, K);   % 扣除参考直流项后的符号能量

for k = 1:K
    for n = 1:Nsym
        idx1 = sym_center(n) - half_win;
        idx2 = idx1 + win_len - 1;

        % 边界保护
        idx1 = max(idx1, 1);
        idx2 = min(idx2, Nt);

        E_sym(n,k) = mean(z_lp(idx1:idx2,k));
        E_sym_remove_ref(n,k) = mean(z_lp_remove_ref(idx1:idx2,k));
    end
end

%% ================= RIS单元平均能量输出 =================
Z = mean(z_lp, 1).';                  % K × 1，未扣除参考项
Z_remove_ref = mean(z_lp_remove_ref, 1).';  % K × 1，扣除参考直流项

%% ================= 可视化：某个RIS单元的时域能量包络 =================
k_plot = 1;

figure;
plot(t*1e6, z_lp(:,k_plot), 'LineWidth', 1.2);
xlabel('Time (\mus)');
ylabel('Detected energy');
title(['Low-pass output after square-law detection, RIS element ', num2str(k_plot)]);
grid on;

%% ================= 可视化：符号化能量检测结果 =================
figure;
stem(1:Nsym, E_sym(:,k_plot), 'filled');
xlabel('Symbol index');
ylabel('Symbol-level energy');
title(['Symbol-level detected energy, RIS element ', num2str(k_plot)]);
grid on;

figure;
stem(1:Nsym, E_sym_remove_ref(:,k_plot), 'filled');
xlabel('Symbol index');
ylabel('Energy after removing reference DC');
title(['Symbol-level energy after reference DC removal, RIS element ', num2str(k_plot)]);
grid on;

%% ================= 可视化：不同RIS单元平均能量 =================
figure;
stem(1:K, Z, 'filled');
xlabel('RIS element index');
ylabel('Average detected energy');
title('Average detected energy across RIS elements');
grid on;