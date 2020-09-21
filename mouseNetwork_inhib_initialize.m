function [simdata,DirPart] = mouseNetwork_inhib_initialize(varies,xrNetcon,Cnoise_coloc,ICdirPath,...
    spks_file,dataCh,plot_rasters,folder,detail,subz,restricts)

datetime = datestr(now,'yyyymmdd-HHMMSS');

ICstruc = dir([ICdirPath '*.mat']);
load(spks_file,'Spks_clean','Spks_masked');

% set(0, 'DefaultFigureVisible', 'off');

h = figure('Position',[50,50,850,690],'visible','off');

% clean cortical noise is stored in varies; save for co-located cases
Cnoise_clean = varies(2).range;

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
    study_dir = fullfile(pwd, folder, [datetime detail], spatialConfig{1});
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
        clean_spks = squeeze(Spks_clean{dataCh}(:,str2double(spatialConfig{1}(2)),:));
        if strcmp(spatialConfig{1}(4),'0')
            data_spks = clean_spks;
        else
            data_spks = squeeze(Spks_masked{dataCh}(:,str2double(spatialConfig{1}(2)),str2double(spatialConfig{1}(4)),:));
        end
    end
    
    % change noise parameter if mixed simulations
    if strcmp(spatialConfig{1}(4),'0') % clean
        varies(2).range = Cnoise_clean;
    else % mixed
        varies(2).range = Cnoise_clean + Cnoise_coloc;
    end

    temp = mouse_network_inhib(study_dir,time_end,varies,xrNetcon,restricts);
   
    [simdata(z).perf,simdata(z).fr,simdata(z).annot,simdata(z).PSTH] = postProcessData(temp,varies,...
        ICdirPath,time_end,study_dir,data_spks,xrNetcon,plot_rasters);
    
    simdata(z).name = ICstruc(z).name;
    
end

Dirparts = strsplit(study_dir,filesep);
DirPart = [filesep fullfile(Dirparts{1:end-1})];

save([DirPart filesep 'summary_results.mat'],'simdata','varies','DirPart')
% close(h);

end