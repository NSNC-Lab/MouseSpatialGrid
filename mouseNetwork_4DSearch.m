clearvars;
close all;

addpath('mechs')
addpath('dependencies')
addpath('eval_scripts')
addpath('genlib')
addpath(genpath('dynasim'))

%%%%%%%% start of user inputs

plot_distances = 0;  % plot VR distances
plot_rasters = 0;   % plot rasters

data_spks_file = '03_30_18_0dB_cleaned(-1,4).mat';
data_perf_file = '03_30_18_0dB_performance.mat';
dataCh = 31;

%%%%%%%% end of user inputs

subject = [extractBefore(data_perf_file,'_performance') '-Ch' num2str(dataCh)];
folder = 'Data-fitting';

%% load experimental data to optimize model to

nCells = 4; %synchronise this variable with mouse_network

load(data_spks_file,'Spks_clean','Spks_masked');
load(data_perf_file,'Max','max_masked','opt_tau','opt_tau_masked');
data_perf = [Max{dataCh}(:);max_masked{dataCh}(:)];

% use tau of best clean location
[~,best_loc] = max(Max{dataCh});
data_tau = opt_tau{dataCh}(best_loc);

%% calculate cortical noise based on spontaneous FR

% data_FR arranged ipsi->contra
for s = 1:4
    data_spks = squeeze(Spks_clean{dataCh}(:,5-s,:));
    FR_r0(s) = mean(cellfun(@(x) sum(x < 0 | x >= 3),data_spks),'all')/2;
    data_FR(s) = mean(cellfun(@(x) sum(x >= 0 & x < 3),data_spks),'all')/3;
end
load('Cnoise_vs_FR0.mat','fit');
Cnoise = (mean(FR_r0)-fit(1))/fit(2);

%% calculate best STRF gain data based on best channel

[~,best_loc] = max(data_perf(1:4));
best_loc = 5-best_loc;    % azimuth indexes are flipped b/w model and data

% extract all STRF gains
gains = dir(fullfile(pwd,'MiceSpatialGrids/ICStim/Mouse/full_grids/BW_0.009 BTM_3.8 t0_0.1 phase0.4985'));
gains(~contains({gains.name},'s30')) = [];

STRF_FR = zeros(length(gains)-1,4);
for g = 1:length(gains)-1
    ICdir = [gains(g).folder filesep gains(g).name];
    ICdirPath = [ICdir filesep];
    ICstruc = dir([ICdirPath '*.mat']);
    
    subz = find(contains({ICstruc.name},'m0.mat')); % sXm0 (target only) cases
    i = 1;
    for z = subz
        % restructure IC spikes
        load([ICdirPath ICstruc(z).name],'avgSpkRate');
        STRF_FR(g,i) = avgSpkRate(best_loc);
        i = i + 1;
    end
end
% find best gain closest to average clean trial FR
[~,ind] = min(abs(mean(STRF_FR,2) - mean(data_FR)));

ICdir = [gains(ind).folder filesep gains(ind).name];
ICdirPath = [ICdir filesep];
ICstruc = dir([ICdirPath '*.mat']);

% back of envelope calculation of upper limit of sum of synaptic
% conductances
gsyn_sum = 0.21*(mean(data_FR) - mean(FR_r0))./(mean(STRF_FR(ind,:)));

%% 4D search
gsyn_range = 0.03; %[0 gsyn_sum/4-0.02:0.01:gsyn_sum/4+0.02];
    
for c = 1:4
   ranges{c} = gsyn_range; 
end

%% varied parameters

varies = struct;

varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;

varies(end+1).conxn = 'C';
varies(end).param = 'noise';
varies(end).range = Cnoise;

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN1';
varies(end).range = ranges{1};

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN2';
varies(end).range = ranges{2};

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN3';
varies(end).range = ranges{3};

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN4';
varies(end).range = ranges{4};

nvaried = {varies(2:end).range};
nvaried = prod(cellfun(@length,nvaried));

[data,DirPart] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,...
    Spks_masked,dataCh,data_tau,plot_distances,plot_rasters,folder,subject,'-4D-search');

%% performance grids

temp = {data.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
%maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
%mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
perf = []; model_FR = [];
for i = 1:length(targetIdx)
    perf(:,5-i) = data(targetIdx(i)).perf.C;
    model_FR(:,5-i) = data(targetIdx(i)).fr.C;
end

[~,MSE_clean] = calcModelPerf(perf,data_perf(1:4)');

loss = MSE_clean(:,1);
frdiffs = abs(mean(model_FR,2) - mean(data_FR(data_FR ~= 0)));
within_thresh = frdiffs <= 5;

% no wts are within FR thresh, just get the absolute minimum
if sum(within_thresh) == 0
   within_thresh = frdiffs == min(frdiffs);
end

makeParallelPlot(data,within_thresh,loss)

[minloss,i] = min(loss(within_thresh));
temp = find(within_thresh);
best_iteration = temp(i);

% make grid of best iteration

makeGrids_bestIteration(data,varies,DirPart,nvaried,data_perf,data_FR,best_iteration);
