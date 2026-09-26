%% Direct validation of the pairwise-error expression over spatial phase maps
clc; clear; close all;
script_dir=fileparts(mfilename('fullpath')); result_dir=fullfile(script_dir,'results');
rng(20261013,'twister'); Ny=8; Nz=8; N=64;
qpsk=exp(1j*(pi/4+(0:3)*pi/2)).'; S=symbol_codebook(qpsk,2);
H=channel_matrix(Ny,Nz,[15;35],[8;-10],[exp(1j*.35);.85*exp(-1j*.55)]);
pool_count=500; code_count=32; phases=zeros(N,pool_count); score=zeros(pool_count,1);
for c=1:pool_count
    if c<=4
        states=(c-1)*ones(N,1);
    elseif c<=20
        [yy,zz]=meshgrid(0:Nz-1,0:Ny-1);
        states=mod((c-4)*yy(:)+(c-3)*zz(:),4);
    else
        dominant=randi(4)-1; dominance=.25+.70*rand;
        states=randi([0,3],N,1); states(rand(N,1)<dominance)=dominant;
    end
    phases(:,c)=(pi/2)*states;
    Mu=abs(H*S+1.5*exp(1j*phases(:,c))).^2;
    scale=sqrt(mean(Mu(:).^2)); D=distance_matrix(Mu);
    score(c)=min(D(D>0))/scale;
end
[~,order]=sort(score); selected=order(round(linspace(1,pool_count,code_count)));
dmin=zeros(code_count,1); pep_mc=zeros(code_count,1); pep_theory=zeros(code_count,1);
snr_dB=0; trials=30000;
for c=1:code_count
    Mu=abs(H*S+1.5*exp(1j*phases(:,selected(c)))).^2;
    scale=sqrt(mean(Mu(:).^2)); sigma=scale/sqrt(10^(snr_dB/10));
    D=distance_matrix(Mu); [d,linear]=min(D(:)); %#ok<ASGLU>
    [i,j]=ind2sub(size(D),linear); dmin(c)=D(i,j)/scale;
    errors=0;
    for t=1:trials
        if rand<.5, truth=i; rival=j; else, truth=j; rival=i; end
        z=Mu(:,truth)+sigma*randn(N,1);
        if norm(z-Mu(:,rival))^2 < norm(z-Mu(:,truth))^2, errors=errors+1; end
    end
    pep_mc(c)=errors/trials;
    pep_theory(c)=.5*erfc(D(i,j)/(2*sigma*sqrt(2)));
end
[dmin,idx]=sort(dmin); pep_mc=pep_mc(idx); pep_theory=pep_theory(idx);
fig=figure('Color','w','Position',[100 100 540 390]);
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman', ...
    'defaultAxesFontSize',11,'defaultLineLineWidth',1.5,'defaultLineMarkerSize',6);
semilogy(dmin,pep_mc,'o','Color',[0 .447 .741],'MarkerFaceColor','w'); hold on;
semilogy(dmin,pep_theory,'k-'); grid on; box on;
xlabel('Normalized minimum energy distance'); ylabel('Worst-pair error probability');
legend('Monte Carlo','Q(d_{min}^{(E)}/(2\sigma_w))','Location','southwest');
exportgraphics(fig,fullfile(result_dir,'fig_spatial_pairwise_error_vs_dmin.eps'), ...
    'ContentType','vector');
exportgraphics(fig,fullfile(result_dir,'fig_spatial_pairwise_error_vs_dmin.png'), ...
    'Resolution',300); savefig(fig,fullfile(result_dir,'fig_spatial_pairwise_error_vs_dmin.fig'));
close(fig);
writetable(table((1:code_count).',dmin,pep_mc,pep_theory,'VariableNames', ...
    {'Code_index','Normalized_dmin','MonteCarlo_PEP','Theoretical_PEP'}), ...
    fullfile(result_dir,'spatial_pairwise_distance_summary.csv'));

function H=channel_matrix(Ny,Nz,td,pd,alpha)
    [yy,zz]=meshgrid(0:Nz-1,0:Ny-1); y=yy(:)-(Nz-1)/2; z=zz(:)-(Ny-1)/2;
    H=zeros(Ny*Nz,numel(td));
    for l=1:numel(td), th=td(l)*pi/180; ph=pd(l)*pi/180;
        H(:,l)=alpha(l)*exp(1j*pi*sin(th).*(y*cos(ph)+z*sin(ph))); end
end
function S=symbol_codebook(a,L)
    c=cell(1,L); [c{:}]=ndgrid(a); S=zeros(L,numel(c{1}));
    for l=1:L, S(l,:)=c{l}(:).'; end
end
function D=distance_matrix(Mu)
    M=size(Mu,2); D=inf(M);
    for i=1:M, for j=i+1:M, D(i,j)=norm(Mu(:,i)-Mu(:,j)); D(j,i)=D(i,j); end, end
end
