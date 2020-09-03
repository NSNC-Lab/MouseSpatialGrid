function [simdata,DirPart] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,...
    Spks_masked,dataCh,plot_rasters,folder,subject,detail,subz,restricts,Cnoise2)

datetime = datestr(now,'yyyymmdd-HHMMSS');

% set(0, 'DefaultFigureVisible', 'off');

h = figure('Position',[50,50,850,690],'visible','off');

Cnoise = varies(2).range;

for z = subz
    % restructure IC spikes
    load([ICdirPath ICstruc(z).name],'t_spiketimes');
    tmat = cellfun(@max,t_spiketimes,'UniformOutput',false);
    tmax = max([tmat{:}]);
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
    
    % change noise parameter if mixed simulations
    if strcmp(spatialConfig{1}(4),'0')
        varies(2).range = Cnoise + Cnoise2;
    else
        varies(2).range = Cnoise;
    end
    

    temp = mouse_network(study_dir,time_end,varies,restricts);
    [simdata(z).perf,simdata(z).fr,simdata(z).annot] = postProcessData(temp,varies,time_end,study_dir,data_spks,plot_rasters);
    
    simdata(z).name = ICstruc(z).name;
    
end

Dirparts = strsplit(study_dir, filesep);
DirPart = [filesep fullfile(Dirparts{1:end-1})];

save([filesep DirPart filesep 'summary_results.mat'],'simdata','varies','DirPart')
% close(h);

end