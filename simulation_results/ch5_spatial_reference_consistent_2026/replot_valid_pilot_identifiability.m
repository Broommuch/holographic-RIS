%% Replot the hardware-consistent columns of the 10,000-trial pilot study
% A fixed reference phase at a representative receiving unit is the scalar
% restriction of the fixed spatial code. Time-varying random/four-phase columns
% from the source study are intentionally excluded.
clc; clear; close all;
script_dir=fileparts(mfilename('fullpath')); result_dir=fullfile(script_dir,'results');
T=readtable(fullfile(result_dir,'pilot_length_source_summary.csv'));
names={'No reference','Constant'}; colors=[.35 .35 .35;0 .447 .741];
fig=figure('Color','w','Position',[100 100 900 340]);
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman', ...
    'defaultAxesFontSize',11,'defaultLineLineWidth',1.5,'defaultLineMarkerSize',6);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for panel=1:2
    nexttile; hold on;
    for k=1:2
        idx=strcmp(T.ReferenceDesign,names{k}); Tk=T(idx,:);
        if panel==1
            plot(Tk.PilotLength,Tk.NoiselessRecoveryProbability,'o-', ...
                'Color',colors(k,:)); ylabel('Recovery probability'); ylim([0 1.03]);
            title('(a) Noiseless multistart recovery');
        else
            plot(Tk.PilotLength,Tk.ChannelNMSEdB,'o-', ...
                'Color',colors(k,:)); ylabel('Channel NMSE (dB)');
            title('(b) Noisy recovery at 12 dB');
        end
    end
    xline(2,':','HandleVisibility','off');
    xline(4,':','HandleVisibility','off');
    xline(6,':','HandleVisibility','off');
    grid on; box on; xlabel('Pilot symbols, T_p'); xlim([1 48]);
end
legend(names,'Location','best','FontSize',9);
exportgraphics(fig,fullfile(result_dir,'fig_spatial_pilot_identifiability.eps'), ...
    'ContentType','vector');
exportgraphics(fig,fullfile(result_dir,'fig_spatial_pilot_identifiability.png'), ...
    'Resolution',300);
savefig(fig,fullfile(result_dir,'fig_spatial_pilot_identifiability.fig')); close(fig);
writetable(T(strcmp(T.ReferenceDesign,'No reference') | ...
    strcmp(T.ReferenceDesign,'Constant'),:), ...
    fullfile(result_dir,'spatial_pilot_identifiability_summary.csv'));
