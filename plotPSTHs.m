function [PSTH,annot] = plotPSTHs(simdata,study_dir,varies)

% figure('position',[0 0 1440 940]);
    
jump = length(find([simdata.IC_IC_trial] == 1));

numTrials = 40;

time_end = length(simdata(1).time);

for vv = 1:jump % for each varied parameter
    
    % subData: all trials for each parameter set
    subData = simdata(vv:jump:length(simdata));
    G_inc(vv) = simdata(vv).C_G_inc;
    tau_ad(vv) = simdata(vv).C_tau_ad;
        
    %% visualize spikes
    Cspks = zeros(numTrials,time_end);
    for i = 1:numTrials
        Cspks(i,:) = subData(i).C_V_spikes;
    end
    
    [PSTH(vv,:)] = makePSTH(Cspks);

    annot(vv,:) = createAnnotStr(simdata(vv));
end
    
% parts = strsplit(study_dir, filesep);
% DirPart = fullfile(parts{1:end-1});
% 
% annotation('textbox',[.675 .85 .2 .1],...
%     'string',annot(vv,:),...
%     'FitBoxToText','on',...
%     'LineStyle','none')
% 
% saveas(gca,[filesep DirPart filesep parts{end} '_PSTH_' num2str(vv) '.tiff'])

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