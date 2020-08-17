function [data,DirPart] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,...
    Spks_masked,dataCh,plot_rasters,folder,subject,detail,best_loc,adaptFlag,allFlag,restrict_vary_flag)

datetime = datestr(now,'yyyymmdd-HHMMSS');

set(0, 'DefaultFigureVisible', 'off');

h = figure('Position',[50,50,850,690]);

if allFlag
    subz = find(~contains({ICstruc.name},'s0'));    % all cases except masker-only
elseif adaptFlag
    subz(1) = find(contains({ICstruc.name},['s' num2str(best_loc) 'm' num2str(best_loc)]));
    subz(2) = cat(2,subz,find(strcmp({ICstruc.name},['s' num2str(best_loc) 'm0.mat']))); % sXm0 (target only) cases
else % fit model to both clean and co-located data
    subz = find(cellfun(@(x) strcmp(x(2),x(4)),{ICstruc.name})); % co-located cases
    subz = cat(2,subz,find(contains({ICstruc.name},'m0.mat'))); % sXm0 (target only) cases
end

for z = subz
    % restructure IC spikes
    load([ICdirPath ICstruc(z).name],'t_spiketimes');
    temp = cellfun(@max,t_spiketimes,'UniformOutput',false);
    tmax = max([temp{:}]);
    spks = zeros(20,4,tmax); %I'm storing spikes in a slightly different way...
    for j = 1:size(t_spiketimes,1) %trials [1:10]
        for k = 1:size(t_spiketimes,2) %neurons [(1:4),(1:4)]
            if k < 5 %song 1
                spks(j,k,round(t_spiketimes{j,k})) = 1;
            else
                spks(j+size(t_spiketimes,1),k-4,round(t_spiketimes{j,k})) = 1;
            end
        end
    end
    
    % save spk file
    spatialConfig = strsplit(ICstruc(z).name,'.');
    study_dir = fullfile(pwd, folder, subject, [datetime detail], spatialConfig{1});
    if exist(study_dir, 'dir')
        rmdir(study_dir, 's');
    end
    mkdir(fullfile(study_dir, 'solve'));
    save(fullfile(study_dir, 'solve','IC_spks.mat'),'spks');
    addpath(fullfile(study_dir, 'solve'))
    
    % call network
    h.Name = ICstruc(z).name;
    time_end = size(spks,3);
    
    % load spikes from data
    if ~contains(spatialConfig{1},'s0')
        if strcmp(spatialConfig{1}(4),'0')
            data_spks = squeeze(Spks_clean{dataCh}(:,str2double(spatialConfig{1}(2)),:));
        else
            data_spks = squeeze(Spks_masked{dataCh}(:,str2double(spatialConfig{1}(2)),str2double(spatialConfig{1}(4)),:));
        end
    else
        data_spks = [];
    end
    
    if adaptFlag
        temp = mouse_network(study_dir,time_end,varies,restrict_vary_flag);
        
        [data(z).PSTH,data(z).tau_ad,data(z).G_inc] = plotPSTHs(temp,study_dir,ICstruc(z).name);
    else % run clean grid or full spatial grid
        temp = mouse_network(study_dir,time_end,varies,restrict_vary_flag);
        
        [data(z).perf,data(z).fr,data(z).annot] = postProcessData(temp,varies,time_end,study_dir,data_spks,plot_rasters);
    end
    
    data(z).name = ICstruc(z).name;
    
end

Dirparts = strsplit(study_dir, filesep);
DirPart = fullfile(Dirparts{1:end-1});

save([filesep DirPart filesep 'summary_results.mat'],'data','varies')
close(h);

end