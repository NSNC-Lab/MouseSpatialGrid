 function [perf,fr,annot,PSTH] = postProcessData(data,varies,ICdirPath,time_end,study_dir,data_spks,plot_rasters)

numTrials = length(varies(1).range);

jump = length(find([data.IC_IC_trial]==1));

STRFgain = extractBetween(ICdirPath,'gain','_2020');

locs = {'90°','45°','0°','-90°'};

tloc = locs(str2double(study_dir(end-2)));
if ~strcmp(study_dir(end),'0') % with masker
    mloc = locs(str2double(study_dir(end)));
    locstr = ['Target @ ' tloc ', Masker @ ' mloc];
else % w/o masker
    locstr = ['Target @ ' tloc];
end


for vv = 1:jump % for each varied parameter
    subData = data(vv:jump:length(data));

    %% visualize spikes
    ICspks = zeros(numTrials,4,time_end);
    %Sspks = zeros(numTrials,4,time_end);
    Rspks = zeros(numTrials,4,time_end);
    Cspks = zeros(numTrials,time_end);
    for i = 1:numTrials
        for j = 1:4
            ICspks(i,j,:) = subData(i).IC_V_spikes(:,j);
            %Sspks(i,j,:) = subData(i).S_V_spikes(:,j);
            Rspks(i,j,:) = subData(i).R_V_spikes(:,j);
        end
        Cspks(i,:) = subData(i).C_V_spikes;
    end

    % plot
    clf;
%     for i = 1:4 %for each spatially directed neuron
%         if i < 4
%             ip = i;
%         else
%             ip = i + 1;
%         end
        
%         if plot_rasters, subplot(4,5,15+ip); end
%         thisRaster = squeeze(ICspks(:,i,:));
%         [perf.IC(i,vv),fr.IC(i,vv)] = calcPCandPlot(thisRaster,time_end,1,plot_rasters,numTrials);        
% 
%         if i==4, ylabel('IC'); end
%         ax = get(gca,'position'); 
%         annotation('textbox',[ax(1)+ax(3)/2 ax(2)-0.02 0 0],...
%             'string',locs{i},...
%             'HorizontalAlignment','center',...
%             'LineStyle','none')
        
%         if plot_rasters, subplot(4,5,10+ip); end
% 
%         thisRaster = squeeze(Sspks(:,i,:));
%         calcPCandPlot(thisRaster,time_end,0,plot_rasters,numTrials);        
%         if i==4, ylabel('S'); end
%         xticklabels([])

%         if plot_rasters, subplot(4,5,5+ip); end
%         thisRaster = squeeze(Rspks(:,i,:));
%         [perf.R(i,vv),fr.R(i,vv)] = calcPCandPlot(thisRaster,time_end,1,plot_rasters,numTrials);
% 
%         if i==4, ylabel('R'); end
%         xticklabels([])
%     end
    
    if plot_rasters %subplot(4,5,3); 
        subplot(2,3,2);
        ylabel('C spikes'); xticklabels([]); 
    end
    [perf.C(vv),fr.C(vv)] = calcPCandPlot(Cspks,time_end,1,plot_rasters,numTrials);     
    
    if plot_rasters, subplot(2,3,1); xticklabels([]); ylabel('Data spikes');
        tempspks = zeros(numel(data_spks),time_end);
        tempspks2 = zeros(numel(data_spks),5000);
        for tt = 1:numel(data_spks)
            temp = round(data_spks{tt}*1000);   % convert spike times from sec to msec
            temp2 = temp; % temp2 includes spontaneous spiking
            
            temp(temp < 0 | temp >= time_end) = [];
            tempspks(tt,temp+1) = 1;    % add ones at spike times
            
            tempspks2(tt,temp2+1001) = 1;    % add ones at spike times for plooting spont.
        end
        
    calcPCandPlot(tempspks,time_end,1,plot_rasters,size(tempspks,1)); 
    
    subplot(2,3,4);
    [PSTH,t_vec] = makePSTH(tempspks);
    plot(t_vec,50*PSTH/size(tempspks,1)); ylim([0 150]);
    xlabel('Time (s)');
    
    end
            
    subplot(2,3,5);
    [PSTH,t_vec] = makePSTH(Cspks);
    plot(t_vec,50*PSTH/40); ylim([0 150]);
    xlabel('Time (s)');
    
    % figure annotations
    annot(vv,:) = createAnnotStr(data(vv),STRFgain);
    
    parts = strsplit(study_dir, filesep);
    DirPart = fullfile(parts{1:end-1});
        
    if plot_rasters
        annotation('textbox',[.675 .85 .2 .1],...
            'string',annot(vv,:),...
            'FitBoxToText','on',...
            'LineStyle','none')
        
        annotation('textbox',[.675 .7 .2 .1],...
            'string', locstr,...
            'FitBoxToText','on',...
            'LineStyle','none')
        
        saveas(gca,[filesep DirPart filesep parts{end} '_v2_' num2str(vv) '.tiff'])
    end
    
end

end

function [pc,fr] = calcPCandPlot(raster,time_end,calcPC,plot_rasters,numTrials)

PCstr = '';

if calcPC
    % spks to spiketimes in a cell array of 20x2
    tau = linspace(1,30,100);
    spkTime = cell(numTrials,1);
    for ii = 1:numTrials, spkTime{ii} = find(raster(ii,:)); end
    spkTime = reshape(spkTime,numTrials/2,2);
    % calculate distance matrix & performance
    distMat = calcvr(spkTime, tau);
    [performance, ~] = calcpc(distMat, numTrials/2, 2, 1,[], 'new');
    pc = mean(max(performance));
    PCstr = ['PC = ' num2str(pc)];
end

if size(raster,2) > 3000
fr = 1000*mean(sum(raster(:,[1:1000,4000:end]),2))/2000;
else
    fr = 1000*mean(sum(raster,2))/3000;
end
%plot
if plot_rasters
    plotSpikeRasterFs(flipud(logical(raster)), 'PlotType','vertline');
    title({PCstr,['FR = ' num2str(fr)]});
    xlim([0 time_end])
    line([0,time_end],[numTrials/2 + 0.5,numTrials/2 + 0.5],'color',[0.3 0.3 0.3])
end
    
end

function annot = createAnnotStr(data,STRFgain)

paramstr = {data(1).varied{2:end}};
gSYNs = []; gs = 1;
i = 1;
for aa = 1:length(paramstr)
    if contains(['data.' paramstr{aa}],'R_C_gSYN')
        gSYNs = cat(2,gSYNs,eval(['data.' paramstr{aa}]));
        gs = gs + 1;
    elseif contains(['data.' paramstr{aa}],'C_noise')
        annot{:,i} = sprintf('%s = %.3f',paramstr{aa},...
            eval(['data.' paramstr{aa}]));
        i = i + 1;
    else
        annot{:,i} = sprintf('%s = %.3f',paramstr{aa},...
            eval(['data.' paramstr{aa}]));
        i = i + 1;
    end
end
% round-up gSYNs so that annotation strings don't have very long numbers
annot{:,end+1} = ['RC_{gSYN} = ' mat2str(round(10000*gSYNs)/10000)];
annot{:,end+1} = ['STRF gain = ' STRFgain{1}];

end