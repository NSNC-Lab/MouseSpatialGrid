function f = mouseNetwork_objective(w)

close all

addpath('mechs')
addpath('dependencies')
addpath('eval_scripts')
addpath('genlib')
addpath(genpath('dynasim'))

plot_grids = 0;  % plot spatial grids
plot_distances = 0;  % plot VR distances
plot_rasters = 0;   % plot rasters

%% define parameters for optimization

if size(w,1) ~= 4
   w = w'; 
end
rcNetcon = w;

%inputChans = find(weights ~= 0);

%% load expeimental data to optimize model to

dataCh = 25;

load('9-21-2016_0dB_removed_trialscleaned(-1,4).mat','Spks_clean','Spks_masked');
load('9-21-2016_0dB_removed_trials_performance.mat','Max','max_masked','opt_tau','opt_tau_masked');
data_perf = [Max{dataCh}(:);max_masked{dataCh}(:)];

% use tau of best clean location
[~,best_loc] = max(Max{dataCh});
data_tau = opt_tau{dataCh}(best_loc);

%% load STRF spikes

researchDrive = 'MiceSpatialGrids/';
ICdir = [researchDrive 'ICStim/Mouse/full_grids//BW_0.009 BTM_3.8 t0_0.1 phase0.4985//s30_STRFgain3.00_20200624-160717'];
ICdirPath = [ICdir filesep];
ICstruc = dir([ICdirPath '*.mat']);
if isempty(ICstruc), error('empty data directory'); end

%% Initialize variables

datetime = datestr(now,'yyyymmdd-HHMMSS');

set(0, 'DefaultFigureVisible', 'off');
h = figure('Position',[50,50,850,690]);
    
%% varied parameters
varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;

varies(end+1).conxn = 'C';
varies(end).param = 'noise';
varies(end).range = 1.6;

nvaried = {varies(2:end).range};
nvaried = prod(cellfun(@length,nvaried));

%% netcons
nCells = 4; %synchronise this variable with mouse_network

% % -90, 0, 45, 90º
% % x-channel inhibition
% irNetcon = zeros(nCells);
% 
% srNetcon = zeros(nCells);

% netCons.irNetcon = irNetcon;
% netCons.srNetcon = srNetcon;
netCons.rcNetcon = rcNetcon;

subz = find(contains({ICstruc.name},'m0.mat')); % sXm0 (target only) cases
%subz = 5:length({ICstruc.name});    % all cases except masker-only
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
    study_dir = fullfile(pwd, 'run', 'optimization', datetime, spatialConfig{1});
    if exist(study_dir, 'dir')
      rmdir(study_dir, 's');
    end
    mkdir(fullfile(study_dir, 'solve'));
    save(fullfile(study_dir, 'solve','IC_spks.mat'),'spks');

    % call network
    h.Name = ICstruc(z).name;
    time_end = size(spks,3);
    
    % load spikes from data
    if ~contains(spatialConfig{1},'s0')
        if strcmp(spatialConfig{1}(4),'0')
            data_spks = squeeze(Spks_clean{dataCh}(:,5-str2double(spatialConfig{1}(2)),:));
        else
            data_spks = squeeze(Spks_masked{dataCh}(:,5-str2double(spatialConfig{1}(2)),5-str2double(spatialConfig{1}(4)),:));
        end
        data_FR(z) = mean(mean(cellfun(@(x) sum(x >= 0 & x < time_end/1000),data_spks)))/time_end*1000;
    else
        data_spks = [];
    end
    
    [data(z).perf, data(z).fr, data(z).annot,~, data(z).VR] = ...
        mouse_network(study_dir,time_end,varies,netCons,plot_rasters,...
        plot_distances,data_spks,data_tau);

    data(z).name = ICstruc(z).name;
end
save([pwd filesep 'run' filesep 'optimization' filesep...
    datetime filesep 'summary_results.mat'],'data')
close(h);

%% performance grids

if plot_grids
    makeGrids(data,varies,DirPart,rcNetcon,nvaried,data_perf,data_FR);
end

temp = {data.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
% maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
% mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));

for i = 1:length(targetIdx)
    perf(i) = data(targetIdx(i)).perf.C;
    fr(i) = data(targetIdx(i)).fr.C;
end
perf = fliplr(perf);

[~,MSE_clean] = calcModelPerf(perf',data_perf(1:4));

f = MSE_clean(1);

end
