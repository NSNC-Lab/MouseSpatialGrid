function [simdata,DirPart] = mouseNetwork_inhibV2_initialize(varies,xrNetcons,Cnoise_coloc,ICdirPath,...
    spks_file,dataCh,plot_rasters,folder,detail,subz,restricts)

datetime = datestr(now,'yyyymmdd-HHMMSS');

ICstruc = dir([ICdirPath '*.mat']);
load(spks_file,'Spks_clean','Spks_masked');

% set(0, 'DefaultFigureVisible', 'off');

h = figure('Position',[50,50,850,690],'visible','off');

% clean cortical noise is stored in varies; save for co-located cases
row = find(strcmp({varies(:).param},'noise'));

Cnoise_clean = varies(row).range;

for z = subz
    
    [~,~,t_end1] = loadICspks(ICdirPath,ICstruc(z).name,folder,datetime,detail);
    [study_dir,spatialConfig,t_end2] = loadICspks...
        (ICdirPath,ICstruc(z+1).name,folder,datetime,detail);
    
    time_end = min([t_end1 t_end2]);
    
    % call network
    h.Name = ICstruc(z).name;
    
    % load spikes from data
    if ~contains(spatialConfig,'s0')
        clean_spks = squeeze(Spks_clean{dataCh}(:,str2double(spatialConfig(2)),:));
        if strcmp(spatialConfig(4),'0')
            data_spks = clean_spks;
        else
            data_spks = squeeze(Spks_masked{dataCh}(:,str2double(spatialConfig(2)),str2double(spatialConfig(4)),:));
        end
    end
    
    % change noise parameter if mixed simulations
    if strcmp(spatialConfig(4),'0') % clean
        varies(row).range = Cnoise_clean;
    else % mixed
        varies(row).range = Cnoise_clean + Cnoise_coloc;
    end

    temp = mouse_network_inhibV2(study_dir,time_end,varies,xrNetcons,restricts);
   
    [simdata(z).perf,simdata(z).fr,simdata(z).annot,simdata(z).PSTH] = postProcessData(temp,varies,...
        ICdirPath,time_end,study_dir,data_spks,xrNetcons,plot_rasters);
    
    simdata(z).name = ICstruc(z).name;
    
end
simdata = simdata(1:2:end);

Dirparts = strsplit(study_dir,filesep);
DirPart = [filesep fullfile(Dirparts{1:end-1})];

save([DirPart filesep 'summary_results.mat'],'simdata','varies','DirPart')
% close(h);

end