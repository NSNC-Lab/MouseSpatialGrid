clearvars;
close all;

addpath('mechs')
addpath('dependencies')
addpath('eval_scripts')
addpath('genlib')
addpath(genpath('dynasim'))
addpath(genpath('ICSimStim'))

%%%%%%%% start of user inputs

data_spks_file = '03_30_18_0dB_cleaned(-1,4).mat';
data_perf_file = '03_30_18_0dB_performance.mat';
dataCh = 31;

%%%%%%%% end of user inputs

subject = [extractBefore(data_perf_file,'_performance') '-Ch' num2str(dataCh)];
folder = 'Data-fitting';

%% load experimental data to optimize model to

nCells = 4; %synchronise this variable with mouse_network

load(data_spks_file,'Spks_clean','Spks_masked','n_cum_clean','n_cum_masked');
load(data_perf_file,'Max','max_masked','opt_tau','opt_tau_masked');

data_perf = Max{dataCh}(:); % clean performance
for t = 1:4
    data_perf(end+1:end+4) = flipud(max_masked{dataCh}(:,t)); % mixed performance
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
    FR_r0(s+4) = mean(cellfun(@(x) sum(x < 0 | x >= 3),data_spks),'all')/2;
    data_FR_masked(s) = mean(cellfun(@(x) sum(x >= 0 & x < 3),data_spks),'all')/3;
end
load('Cnoise_vs_FR0.mat','fit');
Cnoise = (mean(FR_r0)-fit(1))/fit(2);

%% calculate best STRF gain data based on best channel

% fit STRF based on peak inital FR in data
t_vec = -1:1/50:4;
t_inds = t_vec >= 0 & t_vec < 3;

% calculate STRF gain based on average FR
gain = ((data_FR(best_loc)-mean(FR_r0))/12.4832)^(1/0.7172);

ICdir = InputGaussianSTRF_fitting(gain);

ICdirPath = [ICdir filesep];
ICstruc = dir([ICdirPath '*.mat']);

% clean STRF rates
ICs = find(contains({ICstruc.name},'m0') & ~contains({ICstruc.name},'s0'));
for loc = 1:4   % target location
    load([ICdirPath ICstruc(ICs(loc)).name],'avgSpkRate')
    STRF_FR_clean(loc,:) = avgSpkRate;
end

% collocated STRF rates
ICs = find((cellfun(@(x) x(2),{ICstruc.name}) == cellfun(@(x) x(4),{ICstruc.name})));
for loc = 1:4   % target location
    load([ICdirPath ICstruc(ICs(loc)).name],'avgSpkRate')
    STRF_FR_coloc(loc,:) = avgSpkRate;
end

frac_clean = (data_FR(best_loc) - mean(FR_r0)) / mean(STRF_FR_clean(best_loc,:));
frac_coloc = (data_FR_masked(1+5*(best_loc-1)) - mean(FR_r0)) / mean(STRF_FR_coloc(best_loc,:));

gsyn_sum = max([frac_coloc*0.18 frac_clean*0.18]);

varies = struct;

ranges = cell(1,4);
for c = 1:4
    ranges{c} = 0; 
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

subz = find(contains({ICstruc.name},['s' num2str(best_loc) 'm0'])); % co-located cases

[simdata_noise] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,...
    Spks_masked,dataCh,1,folder,subject,'-noise-only',subz,[],0);

%% Estimate cortical noise for co-located trials

varies = struct;

ranges = cell(1,4);
for c = 1:4
    if c == best_loc
        ranges{c} = 0.18*(data_FR(best_loc)-mean(FR_r0))/(STRF_FR_clean(best_loc));
    else
        ranges{c} = 0;
    end
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

subz = find(contains({ICstruc.name},['s' num2str(best_loc) 'm0'])); % co-located cases

[simdata_STRF] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,...
    Spks_masked,dataCh,1,folder,subject,'-STRF-with-noise',subz,[],0);

varies(2).range = Cnoise + [0:0.2:2];

subz = find(contains({ICstruc.name},['s' num2str(best_loc) 'm' num2str(best_loc)])); % co-located cases

[simdata] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,...
    Spks_masked,dataCh,1,folder,subject,'-collocated-noise',subz,[],0);

%% 4D search

temp = {simdata.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content
colocIdx = find((cellfun(@(x) x(2),temp) == cellfun(@(x) x(4),temp)));

[~,ind] = min(abs(simdata(colocIdx).fr.C - data_FR_masked(1 + 5*(best_loc-1))));
Cnoise2 = varies(2).range(ind);

% back of envelope calculation of upper limit of sum of synaptic
% conductances

gsyn_range = 0:0.02:gsyn_sum;

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

[simdata,DirPart] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,Spks_masked,...
    dataCh,plot_rasters,folder,subject,'-4D-search',subz,restricts,Cnoise2);

%% performance grids for 4D search

temp = {simdata.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
colocIdx = find((cellfun(@(x) x(2),temp) == cellfun(@(x) x(4),temp)));
perf_clean = []; model_FR = [];

% clean performance and FR
for i = 1:length(targetIdx)
    perf_clean(:,i) = simdata(targetIdx(i)).perf.C;
    model_FR(:,i) = simdata(targetIdx(i)).fr.C;
end

perf_masked = []; model_FR_masked = [];
% mixed performance and FR
if ~isempty(colocIdx)
for i = 1:length(colocIdx)
    perf_masked(:,i) = simdata(colocIdx(i)).perf.C;
    model_FR_masked(:,i) = simdata(colocIdx(i)).fr.C;
end
end

[~,MSE_clean_perf] = calcModelLoss(perf_clean,data_perf(1:4)');
% calculate loss for FR as a percentage of model FR
[~,MSE_clean_FR] = calcModelLoss(100*(model_FR)./data_FR,100*data_FR./data_FR);

if ~isempty(colocIdx)
    [~,MSE_masked_perf] = calcModelLoss(perf_masked,data_perf([5 10 15 20])');
    % calculate loss for FR as a percentage of model FR
    [~,MSE_masked_FR] = calcModelLoss(100*(model_FR_masked)./data_FR_masked([1 6 11 16]),100*data_FR_masked([1 6 11 16])./data_FR_masked([1 6 11 16]));
    loss = MSE_clean_perf(:,1) + MSE_clean_FR(:,1) + MSE_masked_perf(:,1) + MSE_masked_FR(:,1);
else
    loss = MSE_clean_perf(:,1) + MSE_clean_FR(:,1);
end

frdiffs = abs(mean(model_FR,2) - mean(data_FR(data_FR ~= 0)));
within_thresh = frdiffs <= 5;

% no wts are within FR thresh, just get the absolute minimum
if sum(within_thresh) == 0
   within_thresh = frdiffs == min(frdiffs);
end

% plotFRvsgSYN(model,within_thresh);

% find sets of parameters that are within threshold
[temp,inds] = sort(loss(within_thresh),'ascend');

if sum(within_thresh) >= 5
    best_iterations = zeros(1,5);
    makeParallelPlot(simdata,within_thresh,loss,DirPart);
else
    best_iterations = zeros(1,sum(within_thresh));
end

for i = 1:length(best_iterations)
    best_iterations(i) = find(loss == temp(i));
end

% make grid of best iterations
makeGrids_bestIteration(simdata,DirPart,data_perf,data_FR,data_FR_masked(1:5:end),best_iterations,loss);

%% rerun dynasim and obtain rasters for best iteration

gsyn_strs = cellfun(@str2num,extractAfter({simdata(targetIdx(1)).annot{:,end}},'RC_{gSYN} = '),'UniformOutput',false);
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
%maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
colocIdx = find((cellfun(@(x) x(2),temp) == cellfun(@(x) x(4),temp)));
perf_clean = []; model_FR = [];

% clean performance and FR
for i = 1:length(targetIdx)
    perf_clean(:,i) = simdata(targetIdx(i)).perf.C;
    model_FR(:,i) = simdata(targetIdx(i)).fr.C;
end

perf_masked = []; model_FR_masked = [];
% mixed performance and FR
if ~isempty(colocIdx)
for i = 1:length(colocIdx)
    perf_masked(:,i) = simdata(colocIdx(i)).perf.C;
    model_FR_masked(:,i) = simdata(colocIdx(i)).fr.C;
end
end

[~,MSE_clean_perf] = calcModelLoss(perf_clean,data_perf(1:4)');
% calculate loss for FR as a percentage of model FR
[~,MSE_clean_FR] = calcModelLoss(100*(model_FR)./data_FR,100*data_FR./data_FR);

if ~isempty(colocIdx)
    [~,MSE_masked_perf] = calcModelLoss(perf_masked,data_perf(5:5:end)');
    % calculate loss for FR as a percentage of model FR
    [~,MSE_masked_FR] = calcModelLoss(100*(model_FR_masked)./data_FR_masked(1:5:end),100*data_FR_masked(1:5:end)./data_FR_masked(1:5:end));
    loss = MSE_clean_perf(:,1) + MSE_clean_FR(:,1) + MSE_masked_perf(:,1) + MSE_masked_FR(:,1);
else
    loss = MSE_clean_perf(:,1) + MSE_clean_FR(:,1);
end

makeGrids(simdata,DirPart,data_perf,data_FR,data_FR_masked(1:5:end),loss);
