%% Perfect-CSI detection with a simultaneous spatial four-phase reference code
clc; clear; close all;
script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, 'results');
if ~exist(result_dir, 'dir'), mkdir(result_dir); end
rng(20261011, 'twister');

Ny = 8; Nz = 8; N = Ny * Nz; L = 2;
theta = [15; 35]; phi = [8; -10];
alpha = [exp(1j * 0.35); 0.85 * exp(-1j * 0.55)];
H = channel_matrix(Ny, Nz, theta, phi, alpha);
qpsk = exp(1j * (pi / 4 + (0:3) * pi / 2)).';
S = symbol_codebook(qpsk, L);
M = size(S, 2);
phase = (pi / 2) * mod((0:N-1).', 4);
b = 1.5 * exp(1j * phase);
Mu = abs(H * S + b).^2;
dmat = pairwise_distances(Mu);
dmin = min(dmat(dmat > 0));

snr_dB = -15:2:15;
trials = 30000;
ser_ml = zeros(size(snr_dB));
ser_lin = zeros(size(snr_dB));
ser_gs = zeros(size(snr_dB));
union_bound = zeros(size(snr_dB));
bp_fixed = 25 * ones(size(snr_dB));
bp_cert = zeros(size(snr_dB));

fprintf('\nSpatial-reference perfect-CSI detection\n');
for k = 1:numel(snr_dB)
    sigma2 = mean(Mu(:).^2) / 10^(snr_dB(k) / 10);
    sigma = sqrt(sigma2);
    union_bound(k) = symbol_union_bound(dmat, S, sigma);
    err_ml = 0; err_lin = 0; err_gs = 0; cert_total = 0;
    for t = 1:trials
        idx = randi(M);
        z = Mu(:, idx) + sigma * randn(N, 1);
        [~, ihat] = min(sum((Mu - z).^2, 1));
        err_ml = err_ml + sum(S(:, ihat) ~= S(:, idx));
        s0 = linearized_estimate(z, H, b, qpsk);
        err_lin = err_lin + sum(s0 ~= S(:, idx));
        [shat, used] = gs_detect(z, H, b, qpsk, s0, 25, dmin);
        err_gs = err_gs + sum(shat ~= S(:, idx));
        cert_total = cert_total + used;
    end
    ser_ml(k) = err_ml / (trials * L);
    ser_lin(k) = err_lin / (trials * L);
    ser_gs(k) = err_gs / (trials * L);
    bp_cert(k) = cert_total / trials;
    fprintf('SNR %5.1f: ML %.4g, GS %.4g, bound %.4g, BP %.2f\n', ...
        snr_dB(k), ser_ml(k), ser_gs(k), union_bound(k), bp_cert(k));
end

floor_value = 0.5 / (trials * L);
fig = publication_figure([100 100 900 340]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
semilogy(snr_dB, max(ser_ml, floor_value), 'o-', 'Color', [0 .447 .741]); hold on;
semilogy(snr_dB, max(ser_lin, floor_value), '^-', 'Color', [.466 .674 .188]);
semilogy(snr_dB, max(ser_gs, floor_value), 's--', 'Color', [.85 .325 .098]);
semilogy(snr_dB, max(union_bound, floor_value), 'k-.');
grid on; box on; xlabel('SNR (dB)'); ylabel('SER');
title('(a) Detection accuracy');
legend('ML', 'Linearized', 'Certified GS', 'Union bound', ...
    'Location', 'southwest', 'FontSize', 9);
nexttile;
plot(snr_dB, bp_fixed, 'k--'); hold on;
plot(snr_dB, bp_cert, 'o-', 'Color', [0 .447 .741]);
grid on; box on; xlabel('SNR (dB)'); ylabel('Average GS back-projections');
title('(b) Certified stopping');
legend('Fixed budget', 'Exact-distance certificate', 'Location', 'southwest', ...
    'FontSize', 9);
export_figure(fig, result_dir, 'fig_spatial_certified_gs');

summary = table(snr_dB(:), ser_ml(:), ser_lin(:), ser_gs(:), ...
    union_bound(:), bp_cert(:), 'VariableNames', ...
    {'SNR_dB','ML_SER','Linearized_SER','GS_SER','Union_bound','Certified_BP'});
writetable(summary, fullfile(result_dir, 'spatial_detection_summary.csv'));

save(fullfile(result_dir, 'spatial_detection_results.mat'), 'snr_dB', ...
    'ser_ml','ser_lin','ser_gs','union_bound','bp_cert','H','b','S');

function H = channel_matrix(Ny, Nz, theta_deg, phi_deg, alpha)
    [yy, zz] = meshgrid(0:Nz-1, 0:Ny-1);
    y = yy(:) - (Nz-1)/2; z = zz(:) - (Ny-1)/2;
    H = zeros(Ny*Nz, numel(theta_deg));
    for l = 1:numel(theta_deg)
        th = theta_deg(l)*pi/180; ph = phi_deg(l)*pi/180;
        a = exp(1j*pi*sin(th).*(y*cos(ph)+z*sin(ph)));
        H(:,l) = alpha(l)*a;
    end
end

function S = symbol_codebook(alphabet, L)
    cells = cell(1,L); [cells{:}] = ndgrid(alphabet);
    S = zeros(L, numel(cells{1}));
    for l=1:L, S(l,:) = cells{l}(:).'; end
end

function d = pairwise_distances(Mu)
    M = size(Mu,2); d = inf(M);
    for i=1:M
        for j=i+1:M
            d(i,j)=norm(Mu(:,i)-Mu(:,j)); d(j,i)=d(i,j);
        end
    end
end

function ub = symbol_union_bound(d, S, sigma)
    [L,M] = size(S); ub = 0;
    for i=1:M
        for j=1:M
            if i~=j
                h = sum(S(:,i)~=S(:,j))/L;
                ub = ub + h * 0.5*erfc(d(i,j)/(2*sigma*sqrt(2)));
            end
        end
    end
    ub = min(ub/M, 1);
end

function s = linearized_estimate(z, G, b, alphabet)
    A = 2*[real(conj(b).*G), -imag(conj(b).*G)];
    x = (A'*A + 1e-6*eye(size(A,2))) \ (A'*(z-abs(b).^2));
    L=size(G,2); sc=x(1:L)+1j*x(L+1:end); s=project(sc,alphabet);
end

function [s, used] = gs_detect(z,G,b,alphabet,s,maxit,dmin)
    used=0;
    if norm(z-abs(G*s+b).^2) < dmin/2, return; end
    sc=s;
    for it=1:maxit
        u=G*sc+b; y=sqrt(max(z,0)).*exp(1j*angle(u));
        sc=(G'*G+1e-6*eye(size(G,2)))\(G'*(y-b));
        s=project(sc,alphabet); used=it;
        if norm(z-abs(G*s+b).^2) < dmin/2, return; end
    end
end

function s=project(x,a)
    s=zeros(size(x));
    for k=1:numel(x), [~,i]=min(abs(x(k)-a)); s(k)=a(i); end
end

function fig=publication_figure(pos)
    fig=figure('Color','w','Position',pos);
    set(groot,'defaultAxesFontName','Times New Roman', ...
        'defaultTextFontName','Times New Roman','defaultAxesFontSize',11, ...
        'defaultLineLineWidth',1.5,'defaultLineMarkerSize',6);
end

function export_figure(fig,dir,name)
    exportgraphics(fig,fullfile(dir,[name '.eps']),'ContentType','vector');
    exportgraphics(fig,fullfile(dir,[name '.png']),'Resolution',300);
    savefig(fig,fullfile(dir,[name '.fig'])); close(fig);
end
