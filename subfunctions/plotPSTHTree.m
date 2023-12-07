function plotPSTHTree(data,configName,options)
% calculate performance and FR for *a single spot on the spatial grid*
% input:
%   data structure with length = #sims, containing voltage information of
%   [time x spatial channel]. Time dimension contains information from each
%   spot in the spatial grid. Extract this information with
%   options.trialStart and options.trialEnd.
% output:
%   data
%       spatial channels,trials,time
% no plotting function for now

trial_length = options.trial_length;
dt = options.dt;

jump = length(data.spks.On);
SpatAttention = options.SpatialAttention;

% network properties
popNames = fieldnames(data.fr);
nPops = numel(popNames); 

nChans = length(fieldnames(data.perf.On));
annotConfig = configName(end-3:end);

% visualize spikes for specified populations
if ~isfield(options,'subPops'), options.subPops = popNames; end

% in same order as popNames
if nPops == 8 % no top down
    subplot_locs = [11 12 7 9 8 4 5 2]; %[11 7 8 4 5 9 6 2];
elseif nPops == 9 % no top down
    subplot_locs = [11 10 7 8 4 5 9 6 2];
elseif nPops == 10 % no C neuron
    subplot_locs = [14 16 10 12 9 11 6 8 5 7 ];
elseif nPops == 14 % 2 X neurons, 1 C neuron, 1 TD neuron
    subplot_locs = [14 16 10 12 9 11 6 8 5 7 1 13 15 2];
else
    subplot_locs = [14 10 11 7 8 12 9 1 2 5];
end

figure('unit','inches','position',[6 3 6 5]);

for vv = 1:jump % for each varied parameter
        
    subSpks = data.spks(vv);
    subPC = data.perf(vv);
    subFR = data.fr(vv);
    
    for ch = 1:nChans

        for currentPop = 1:nPops

            if strcmp(popNames{currentPop},'C'), ch_num = 1;
            else, ch_num = ch; end

            plotSubPSTH(subSpks.(popNames{currentPop}).(['channel' num2str(ch_num)]),...
                subPC.(popNames{currentPop}).(['channel' num2str(ch_num)]),...
                subFR.(popNames{currentPop}).(['channel' num2str(ch_num)]),...
                trial_length,dt,popNames{currentPop},subplot_locs(currentPop),SpatAttention);
        end

        % make ylim the same across all subplots
        psths = findobj(gcf,'type','line');
        ymax = max([psths.YData]);
        plts = findobj(gcf,'type','axes');
        for p = 1:length(plts)
            plts(p).YLim = [0 ymax];
        end

        figName = sprintf('%s_CH%i_set%s_PSTH',configName,ch,num2str(vv));
        
        annotation('textbox',[.6 .82 .1 .1], ...
            'String',[annotConfig, ', CH' num2str(ch) ],'EdgeColor','none','FontSize',20)
        
        saveas(gcf,[figName '.png']);
        % savefig(gcf,[figName '.fig']);
        clf;
    end

end

end

function plotSubPSTH(raster,pc,fr,trial_length,dt,unit,subplot_loc,SpatAttention)

PCstr = ['PC = ' num2str(round(pc)) '%'];

% ind2sub counts down per column first, 
if SpatAttention
    [c,r] = ind2sub([3 5],subplot_loc);
    r = 6-r;
        ypos = 0.06 + 0.2*(r-1);
y = 0.10;

else
    [c,r] = ind2sub([4 4],subplot_loc);
    r = 5-r;
    ypos = 0.06 + 0.215*(r-1);
y = 0.12;

end

xpos = 0.06 + 0.23*(c-1);
x = 0.22;

subplot('position',[xpos ypos x y]); hold on;

% plot target responses on top of each other

t_bin = 20; % in [ms] 
t_vec = (0:t_bin:trial_length)/dt;

colors = {'k','r'};
for tid = 1:2
    [~,spk_inds] = find(raster((1:10) + 10*(tid-1),:));
    PSTH = histcounts(spk_inds,t_vec);
    plot(dt*t_vec(1:end-1),PSTH,colors{tid});
end
xlim([0 t_vec(end)*dt]);
title({unit,[PCstr,[', FR = ' num2str(fr)]]},'fontweight','normal','fontsize',8); set(gca,'xtick',[],'ytick',[])

end
