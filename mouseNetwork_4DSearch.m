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
data_perf = Max{dataCh}(:);
for t = 1:4
    data_perf(end+1:end+4) = flipud(max_masked{dataCh}(:,t));
end
% data_perf(1:4) - clean perf contra->ipsi
% data_perf(5:end) - target(m m m m) target (m m m m) etc.
        % target - contra->ipsi
        % masker - contra->ipsi

% use tau of best clean location
[~,best_loc] = max(Max{dataCh});
data_tau = opt_tau{dataCh}(best_loc);

%% calculate cortical noise based on spontaneous FR

% data_spks arranged contra->ipsi
FR_r0 = zeros(1,4);
data_FR = zeros(1,4);
for s = 1:4
    data_spks = squeeze(Spks_clean{dataCh}(:,s,:));
    FR_r0(s) = mean(cellfun(@(x) sum(x < 0 | x >= 3),data_spks),'all')/2;
    data_FR(s) = mean(cellfun(@(x) sum(x >= 0 & x < 3),data_spks),'all')/3;
end

% data_FR_masked arranged: 4 inds for target at -90, 4 inds at 0, etc.
% similar arrangement as data_FR for clean
for s = 1:16
    ti = ceil(s/4);
    mi = mod(s,4);
    if mi == 0
       mi = 4; 
    end
   data_spks = squeeze(Spks_masked{dataCh}(:,ti,mi,:));
   data_FR_masked(s) = mean(cellfun(@(x) sum(x >= 0 & x < 3),data_spks),'all')/3;
end
load('Cnoise_vs_FR0.mat','fit');
Cnoise = (mean(FR_r0)-fit(1))/fit(2);

%% calculate best STRF gain data based on best channel

% extract all STRF gains
gains = dir(fullfile(pwd,'STRFs/BW_0.009 BTM_3.8 t0_0.1 phase0.4985'));
gains(~contains({gains.name},'s30')) = [];

STRF_FR = zeros(length(gains),4);
for g = 1:length(gains)
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
    
ranges = cell(1,4);
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
varies(end).range = 0;%ranges{1};

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN2';
varies(end).range = 0.03;%ranges{2};

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN3';
varies(end).range = 0.06;%ranges{3};

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN4';
varies(end).range = 0.03;%ranges{4};

nvaried = {varies(2:end).range};
nvaried = prod(cellfun(@length,nvaried));

[data,DirPart] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,...
    Spks_masked,dataCh,data_tau,plot_distances,plot_rasters,folder,subject,'-4D-search');

%% performance grids

temp = {data.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
%maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
perf = []; model_FR = [];

% clean performance and FR
for i = 1:length(targetIdx)
    perf(:,i) = data(targetIdx(i)).perf.C;
    model_FR(:,i) = data(targetIdx(i)).fr.C;
end

perf_masked = []; model_FR_masked = [];
% mixed performance and FR
for i = 1:length(mixedIdx)
    perf_masked(:,i) = data(mixedIdx(i)).perf.C;
    model_FR_masked(:,i) = data(mixedIdx(i)).fr.C;
end

[~,MSE_clean_perf] = calcModelLoss(perf,data_perf(1:4)');
% calculate loss for FR as a percentage of model FR
[~,MSE_clean_FR] = calcModelLoss(100*(model_FR)./data_FR,100*data_FR./data_FR);

[~,MSE_masked_perf] = calcModelLoss(perf_masked,data_perf(5:end)');
% calculate loss for FR as a percentage of model FR
[~,MSE_masked_FR] = calcModelLoss(100*(model_FR_masked)./data_FR_masked,100*data_FR_masked./data_FR_masked);

loss = MSE_clean_perf(:,1) + MSE_clean_FR(:,1) + MSE_masked_perf(:,1) + MSE_masked_FR(:,1);

frdiffs = abs(mean(model_FR,2) - mean(data_FR(data_FR ~= 0)));
within_thresh = frdiffs <= 5;

% no wts are within FR thresh, just get the absolute minimum
if sum(within_thresh) == 0
   within_thresh = frdiffs == min(frdiffs);
end

makeParallelPlot(data,within_thresh,loss)

[temp,inds] = sort(loss(within_thresh),'ascend');

best_iterations = zeros(1,5);
for i = 1:5
    best_iterations(i) = find(loss == temp(i));
end

% make grid of best iterations
makeGrids_bestIteration(data,varies,DirPart,data_perf,data_FR,best_iterations,loss);
