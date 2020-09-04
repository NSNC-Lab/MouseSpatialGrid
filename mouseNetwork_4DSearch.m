% Run this part after mouseNetwork_premanual

% use this code to let Matlab find the best fit for you (takes a long time)

clearvars;
close all;

addpath('mechs')
addpath('dependencies')
addpath('eval_scripts')
addpath('genlib')
addpath(genpath('dynasim'))
addpath(genpath('ICSimStim'))

%%%%%%%% start of user inputs

spks_file = '03_30_18_0dB_cleaned(-1,4).mat';
perf_file = '03_30_18_0dB_performance.mat';
dataCh = 31;

fitFlag = 1;    % if only want to generate clean and co-lcoated spots, =1

%%%%%%%% end of user inputs

subject = [extractBefore(perf_file,'_performance') '-Ch' num2str(dataCh)];
folder = ['Data-fitting' filesep subject];

load(fullfile(folder,'default_parameters.mat'));

ICdirPath = [ICdir filesep];
ICstruc = dir([ICdirPath '*.mat']);

%% varied parameters

varies = struct;

ranges = cell(1,4);
for n = 1:4
    ranges{n} = 0:0.02:gsyn_sum;
end

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

% for now: last varies field has to be restriction
restricts(1).conxn = 'R->C';
restricts(1).param = {'gSYN1','gSYN2','gSYN3','gSYN4'};
restricts(1).range = [gsyn_sum/2 gsyn_sum];

gsynchans = {'gSYN1','gSYN2','gSYN3','gSYN4'};

restricts(end+1).conxn = 'R->C';
restricts(end).param = gsynchans(best_loc);
restricts(end).range = [gsyn_sum/2 gsyn_sum];

plot_rasters = 0;

subz = find(cellfun(@(x) strcmp(x(2),x(4)),{ICstruc.name})); % co-located cases
subz = cat(2,subz,find(contains({ICstruc.name},'m0.mat'))); % sXm0 (target only) cases

[simdata,DirPart] = mouseNetwork_initialize(varies,Cnoise2,ICdirPath,spks_file,dataCh,plot_rasters,...
    folder,'-4D-search',subz,restricts);

%% performance grids for 4D search

nvaried = size(simdata(1).annot,1);

temp = {simdata.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
colocIdx = find((cellfun(@(x) x(2),temp) == cellfun(@(x) x(4),temp)));

perf_clean = zeros(nvaried,4); model_FR = zeros(nvaried,4);
% clean performance and FR
for i = 1:length(targetIdx)
    perf_clean(:,i) = simdata(targetIdx(i)).perf.C;
    model_FR(:,i) = simdata(targetIdx(i)).fr.C;
end

perf_masked = zeros(nvaried,4);
% mixed performance and FR
if ~isempty(colocIdx)
    for i = 1:length(colocIdx)
        perf_masked(:,i) = simdata(colocIdx(i)).perf.C;
    end
end

[~,MSE_clean_perf] = calcModelLoss(perf_clean,data_perf(1:4)');

if ~isempty(colocIdx)
    [~,MSE_masked_perf] = calcModelLoss(perf_masked,data_perf([5 10 15 20])');
    loss = MSE_clean_perf(:,1)+ MSE_masked_perf(:,1);
else
    loss = MSE_clean_perf(:,1);
end

% if no sets are within FR threshold, just find the set with the lowest
% loss

frdiffs = abs(mean(model_FR,2) - mean(data_FR(data_FR ~= 0)));
within_thresh = frdiffs <= 5;

if sum(within_thresh) == 0
   within_thresh = frdiffs == min(frdiffs);
end

% find sets of parameters that are within threshold
[temp,inds] = sort(loss(within_thresh),'ascend');

if sum(within_thresh) >= 5
    best_iterations = zeros(1,5);
    % makeParallelPlot(simdata,within_thresh,loss,DirPart);
else
    best_iterations = zeros(1,sum(within_thresh));
end

for i = 1:length(best_iterations)
    best_iterations(i) = find(loss == temp(i));
end

% make grid of best iterations
makeGrids_fitting(simdata,DirPart,data_perf,data_FR,data_FR_masked(1:5:end));

%% rerun dynasim and obtain rasters for best iteration

gsyn_strs = cellfun(@str2num,extractAfter(simdata(targetIdx(1)).annot(:,end-1),'RC_{gSYN} = '),'UniformOutput',false);
best_gsyns = gsyn_strs{best_iterations(1)};

varies = struct;

varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;

varies(end+1).conxn = 'C';
varies(end).param = 'noise';
varies(end).range = Cnoise;

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN1';
varies(end).range = best_gsyns(1);

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN2';
varies(end).range = best_gsyns(2);

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN3';
varies(end).range = best_gsyns(3);

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN4';
varies(end).range = best_gsyns(4);

subz = find(~contains({ICstruc.name},'s0')); % run all spots on grid
restricts = [];
plot_rasters = 1;

[simdata,DirPart] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,...
    Spks_masked,dataCh,plot_rasters,folder,subject,'-best-iteration',subz,restricts,Cnoise2);

%% make full spatial grid

temp = {simdata.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
colocIdx = find((cellfun(@(x) x(2),temp) == cellfun(@(x) x(4),temp)));

perf_clean = zeros(1,4);
% clean performance and FR
for i = 1:length(targetIdx)
    perf_clean(i) = simdata(targetIdx(i)).perf.C;
end

perf_masked = zeros(1,4);
% mixed performance and FR
if ~isempty(colocIdx)
    for i = 1:length(colocIdx)
        perf_masked(i) = simdata(colocIdx(i)).perf.C;
    end
end

[~,MSE_clean_perf] = calcModelLoss(perf_clean,data_perf(1:4)');

if ~isempty(colocIdx)
    [~,MSE_masked_perf] = calcModelLoss(perf_masked,data_perf(5:5:end)');
    loss = MSE_clean_perf(:,1) + MSE_masked_perf(:,1);
else
    loss = MSE_clean_perf(:,1);
end

makeGrids(simdata,DirPart,data_perf,data_FR,data_FR_masked(1:5:end),loss);
