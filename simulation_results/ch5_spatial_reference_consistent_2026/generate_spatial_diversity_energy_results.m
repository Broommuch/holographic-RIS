%% Spatial diversity, angular separation, and reference-energy balance
clc; clear; close all;
script_dir=fileparts(mfilename('fullpath'));
result_dir=fullfile(script_dir,'results');
if ~exist(result_dir,'dir'), mkdir(result_dir); end
rng(20261012,'twister');
qpsk=exp(1j*(pi/4+(0:3)*pi/2)).'; S=symbol_codebook(qpsk,2);
alpha=[exp(1j*.35); .85*exp(-1j*.55)];

%% Number of spatial observations
Nvec=[4 8 16 32 64 128 256 512 1024];
dminN=zeros(size(Nvec)); serN=zeros(size(Nvec)); snrN=-12; trialsN=10000;
for k=1:numel(Nvec)
    [Ny,Nz]=near_square(Nvec(k));
    H=channel_matrix(Ny,Nz,[15;35],[8;-10],alpha);
    b=1.5*exp(1j*(pi/2)*mod((0:Nvec(k)-1).',4));
    Mu=abs(H*S+b).^2; d=pairwise_distances(Mu);
    scale=sqrt(mean(Mu(:).^2)); dminN(k)=min(d(d>0))/scale;
    sigma=scale/sqrt(10^(snrN/10)); errors=0;
    for t=1:trialsN
        idx=randi(size(S,2)); z=Mu(:,idx)+sigma*randn(Nvec(k),1);
        [~,ih]=min(sum((Mu-z).^2,1)); errors=errors+sum(S(:,ih)~=S(:,idx));
    end
    serN(k)=errors/(2*trialsN);
    fprintf('N %4d: normalized dmin %.3f, SER %.4g\n',Nvec(k),dminN(k),serN(k));
end
fig=publication_figure([100 100 900 340]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; semilogx(Nvec,dminN,'o-','Color',[0 .447 .741]); grid on; box on;
xlabel('Number of receiving units, N_r'); ylabel('Normalized d_{min}^{(E)}');
title('(a) Decision margin');
nexttile; loglog(Nvec,max(serN,.5/(2*trialsN)),'s--','Color',[.85 .325 .098]);
grid on; box on; xlabel('Number of receiving units, N_r'); ylabel('ML SER');
title('(b) Detection reliability');
export_figure(fig,result_dir,'fig_spatial_unit_diversity');
writetable(table(Nvec(:),dminN(:),serN(:),'VariableNames', ...
    {'Num_units','Normalized_dmin','ML_SER'}), ...
    fullfile(result_dir,'spatial_unit_diversity.csv'));

%% Angular separation
sep=[.25 .5 1 2 4 8 12 16 24 32];
sigmin=zeros(size(sep)); serSep=zeros(size(sep)); trialsSep=20000; snrSep=-6;
for k=1:numel(sep)
    H=channel_matrix(8,8,[20;20+sep(k)],[5;5],alpha);
    V=H./reshape(alpha,1,[]); sigmin(k)=min(svd(V))/sqrt(64);
    b=1.5*exp(1j*(pi/2)*mod((0:63).',4)); Mu=abs(H*S+b).^2;
    scale=sqrt(mean(Mu(:).^2)); sigma=scale/sqrt(10^(snrSep/10)); errors=0;
    for t=1:trialsSep
        idx=randi(size(S,2)); z=Mu(:,idx)+sigma*randn(64,1);
        [~,ih]=min(sum((Mu-z).^2,1)); errors=errors+sum(S(:,ih)~=S(:,idx));
    end
    serSep(k)=errors/(2*trialsSep);
end
fig=publication_figure([100 100 900 340]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; semilogx(sep,sigmin,'o-','Color',[0 .447 .741]); grid on; box on;
xlabel('Angular separation (degree)'); ylabel('\sigma_{min}(V)/\surdN_r');
title('(a) Manifold conditioning');
nexttile; semilogx(sep,max(serSep,.5/(2*trialsSep)),'s--','Color',[.85 .325 .098]);
set(gca,'YScale','log'); grid on; box on; xlabel('Angular separation (degree)');
ylabel('Perfect-CSI ML SER'); title('(b) Detection reliability');
export_figure(fig,result_dir,'fig_spatial_angular_separation');
writetable(table(sep(:),sigmin(:),serSep(:),'VariableNames', ...
    {'Separation_deg','Normalized_sigma_min','ML_SER'}), ...
    fullfile(result_dir,'spatial_angular_separation.csv'));

%% Energy balance: ensemble of fixed balanced spatial phase maps
H0=channel_matrix(8,8,[15;35],[8;-10],alpha);
field0=H0*S; field0=field0/sqrt(mean(abs(field0(:)).^2));
rsr_dB=-24:3:24; Ps=1./(1+10.^(rsr_dB/10)); Pb=1-Ps;
num_codes=24; trials_code=1500; balance_snr=[-6 -3 0];
dmean=zeros(size(rsr_dB)); dstd=zeros(size(rsr_dB)); serRatio=zeros(numel(rsr_dB),numel(balance_snr));
for c=1:num_codes
    states=repmat((0:3).',16,1); states=states(randperm(64));
    unitb=exp(1j*(pi/2)*states);
    k0=find(rsr_dB==0,1);
    Mu0=abs(sqrt(Ps(k0))*field0+sqrt(Pb(k0))*unitb).^2;
    balanced_scale=sqrt(mean(Mu0(:).^2));
    for k=1:numel(rsr_dB)
        Mu=abs(sqrt(Ps(k))*field0+sqrt(Pb(k))*unitb).^2;
        d=pairwise_distances(Mu); values(c,k)=min(d(d>0)); %#ok<SAGROW>
        for sidx=1:numel(balance_snr)
            sigma=balanced_scale/sqrt(10^(balance_snr(sidx)/10)); errors=0;
            for t=1:trials_code
                idx=randi(size(S,2)); z=Mu(:,idx)+sigma*randn(64,1);
                [~,ih]=min(sum((Mu-z).^2,1)); errors=errors+sum(S(:,ih)~=S(:,idx));
            end
            serValues(c,k,sidx)=errors/(2*trials_code); %#ok<SAGROW>
        end
    end
end
dmean=mean(values,1); dstd=std(values,0,1); serRatio=squeeze(mean(serValues,1));
fig=publication_figure([100 100 900 340]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; errorbar(rsr_dB,dmean,dstd,'o-','Color',[0 .447 .741], ...
    'CapSize',4); grid on; box on; xlabel('RSR (dB)');
ylabel('Ensemble-mean d_{min}^{(E)}'); title('(a) Spatial-code separation');
nexttile; styles={'o-','s--','d-.'}; colors=[0 .447 .741;.85 .325 .098;.929 .694 .125];
for sidx=1:numel(balance_snr)
    semilogy(rsr_dB,max(serRatio(:,sidx),.5/(2*trials_code*num_codes)), ...
        styles{sidx},'Color',colors(sidx,:)); hold on;
end
grid on; box on; xlabel('RSR (dB)'); ylabel('ML SER');
title('(b) Detection reliability');
legend(compose('Balanced-case SNR = %d dB',balance_snr),'Location','best','FontSize',8);
export_figure(fig,result_dir,'fig_spatial_reference_energy_balance');
writetable(table(rsr_dB(:),dmean(:),dstd(:),serRatio(:,1),serRatio(:,2),serRatio(:,3), ...
    'VariableNames',{'RSR_dB','Mean_dmin','Std_dmin','SER_m6dB','SER_m3dB','SER_0dB'}), ...
    fullfile(result_dir,'spatial_reference_energy_balance.csv'));
save(fullfile(result_dir,'spatial_diversity_energy_results.mat'),'Nvec','dminN','serN', ...
    'sep','sigmin','serSep','rsr_dB','dmean','dstd','serRatio');

function [a,b]=near_square(N)
    a=floor(sqrt(N)); while mod(N,a)~=0, a=a-1; end; b=N/a;
end
function H=channel_matrix(Ny,Nz,td,pd,alpha)
    [yy,zz]=meshgrid(0:Nz-1,0:Ny-1); y=yy(:)-(Nz-1)/2; z=zz(:)-(Ny-1)/2;
    H=zeros(Ny*Nz,numel(td));
    for l=1:numel(td)
        th=td(l)*pi/180; ph=pd(l)*pi/180;
        H(:,l)=alpha(l)*exp(1j*pi*sin(th).*(y*cos(ph)+z*sin(ph)));
    end
end
function S=symbol_codebook(a,L)
    c=cell(1,L); [c{:}]=ndgrid(a); S=zeros(L,numel(c{1}));
    for l=1:L, S(l,:)=c{l}(:).'; end
end
function d=pairwise_distances(Mu)
    M=size(Mu,2); d=inf(M);
    for i=1:M, for j=i+1:M, d(i,j)=norm(Mu(:,i)-Mu(:,j)); d(j,i)=d(i,j); end, end
end
function fig=publication_figure(pos)
    fig=figure('Color','w','Position',pos); set(groot,'defaultAxesFontName','Times New Roman', ...
        'defaultTextFontName','Times New Roman','defaultAxesFontSize',11, ...
        'defaultLineLineWidth',1.5,'defaultLineMarkerSize',6);
end
function export_figure(fig,dir,name)
    exportgraphics(fig,fullfile(dir,[name '.eps']),'ContentType','vector');
    exportgraphics(fig,fullfile(dir,[name '.png']),'Resolution',300);
    savefig(fig,fullfile(dir,[name '.fig'])); close(fig);
end
