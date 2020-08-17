function [PSTH,tau_ad,G_inc] = plotPSTHs(data,study_dir,name)

figure('position',[0 0 1440 940]);
    
jump = length(find([data.IC_IC_trial] == 1));

numTrials = 40;

time_end = length(data(1).time);

for vv = 1:jump % for each varied parameter
    
    % subData: all trials for each parameter set
    subData = data(vv:jump:length(data));
    G_inc(vv) = data(vv).C_G_inc;
    tau_ad(vv) = data(vv).C_tau_ad;
        
    %% visualize spikes
    Cspks = zeros(numTrials,time_end);
    for i = 1:numTrials
        Cspks(i,:) = subData(i).C_V_spikes;
    end
    
    [PSTH(vv,:,1)] = makePSTH(Cspks(1:numTrials/2,:));
    [PSTH(vv,:,2)] = makePSTH(Cspks(numTrials/2+1:end,:));

end
    
% plot PSTHs
PSTH(tau_ad == 0,:,:) = [];

G_inc_plot = unique(G_inc);
tau_ad_plot = unique(tau_ad); tau_ad_plot(tau_ad_plot == 0) = [];

gi = 0;
ti = 0;

for vv = 1:size(PSTH,1)
    
    % G_inc increases down each row, tau_ad increases to the right
    subplot(length(G_inc_plot),length(tau_ad_plot),vv)
    maxFR = max(PSTH(vv,:,:),[],'all');

    plot(1:20:time_end,PSTH(vv,:,1),'r')
    hold on;
    plot([1 time_end],[maxFR maxFR],'k');
    plot(1:20:time_end,PSTH(vv,:,2) + maxFR,'b')
    ylim([0 maxFR*2])
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    
    if mod(vv,length(tau_ad_plot)) == 1
        gi = gi + 1;
        ylabel(num2str(G_inc_plot(gi)));
    end
    
    if vv > (length(G_inc_plot)-1)*length(tau_ad_plot)
        ti = ti + 1;
        xlabel(num2str(tau_ad_plot(ti)));
    end
end

h1 = axes(gcf,'visible','off');
h1.XLabel.Visible = 'on';
h1.YLabel.Visible = 'on';
xlabel(h1,'Adaptation decay');
ylabel(h1,'Adaptation strength');

sgtitle(name);

saveas(gcf,fullfile(study_dir,['adaptation_PSTHs.tiff']));
clf;

end

function [PSTH] = makePSTH(raster)

% calculate PSTH of model results
t_vec = 1:20:size(raster,2);
temp = sum(raster);
for t = 1:length(t_vec)-1
    PSTH(t) = sum(temp(t_vec(t):t_vec(t+1)));
end
PSTH(end+1) = sum(temp(t_vec(end):end));

end