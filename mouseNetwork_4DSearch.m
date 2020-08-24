clearvars;
close all;

addpath('mechs')
addpath('dependencies')
addpath('eval_scripts')
addpath('genlib')
addpath(genpath('dynasim'))
addpath(genpath('ICSimStim'))

%%%%%%%% start of user inputs

data_spks_file = '9-21-2016_0dB_removed_trials_cleaned(-1,4).mat';
data_perf_file = '9-21-2016_0dB_removed_trials_performance.mat';
dataCh = 25;

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

% Multiply PSTHs by scale factors to get FR (multiply average by # bins
% w/in 1 ms (50) and divide by # trials)
scale_model = 50/20;
scale_data = 50/size(Spks_clean{dataCh},1);

%% calculate cortical noise based on spontaneous FR

onset = 0.05;    % used for FR calculations

% data_spks arranged contra->ipsi
FR_r0 = zeros(1,4);
data_FR = zeros(1,4);
for s = 1:4
    data_spks = squeeze(Spks_clean{dataCh}(:,s,:));
    FR_r0(s) = mean(cellfun(@(x) sum(x < 0 | x >= 3),data_spks),'all')/2;
    data_FR(s) = mean(cellfun(@(x) sum(x >= 0 & x < 3),data_spks),'all')/3;
end

data_PSTH = squeeze(n_cum_clean{dataCh}(best_loc,:)');
data_PSTH = cat(1,data_PSTH{:});

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

data_PSTH_masked = squeeze(n_cum_masked{dataCh}(best_loc,best_loc,:));
data_PSTH_masked = cat(1,data_PSTH_masked{:});


%% calculate best STRF gain data based on best channel

% fit STRF based on peak inital FR in data
t_vec = -1:1/50:4;
t_inds = t_vec >= 0 & t_vec < 3;

FR_mixed_data = scale_data*mean(data_PSTH_masked(:,t_inds),1);
data_onset = max(FR_mixed_data(1:onset*1000/50));

% calculate STRF gain based on data onset

% gain = data_onset/33.3441;
% gain = data_onset/21.5190;  % for co-located data
gain = data_onset/32.2443

ICdir = InputGaussianSTRF_fitting(gain);

ICdirPath = [ICdir filesep];
ICstruc = dir([ICdirPath '*.mat']);

z = find(contains({ICstruc.name},['s' num2str(best_loc) 'm' num2str(best_loc)]));

STRF = load(fullfile(ICdirPath,ICstruc(z).name),'t_spiketimes');
t_vec = 1:20:2986;

temp = [];
for t = 1:20
temp = cat(1,temp,[STRF.t_spiketimes{t,best_loc};STRF.t_spiketimes{t,best_loc+4}]);
end

STRF_PSTH = zeros(size(t_vec));
for t = 1:length(t_vec)-1
   STRF_PSTH(t) = sum(temp >= t_vec(t) & temp < t_vec(t+1));
end
STRF_PSTH(end) = sum(temp >= t_vec(end));

max_PSTH = max([STRF_PSTH,FR_mixed_data]);

STRF_FR = mean(STRF_PSTH);

% plot data vs STRF psths - masked location
figure('visible','on');
subplot(2,1,1); plot(0:0.02:0.02*(length(FR_mixed_data)-1),FR_mixed_data)
title('Experimental data, co-located at best location')
line([onset onset],[0 max_PSTH],'color','r');
subplot(2,1,2); plot(0:0.02:0.02*(length(STRF_PSTH)-1),STRF_PSTH);
line([onset onset],[0 max_PSTH],'color','r');
title('STRF, co-located at best location')

saveas(gcf,fullfile(folder,subject,'data_vs_STRF_PSTH.tiff'));
close;

steady_state_FR = mean(scale_data*mean(data_PSTH_masked(:,1.5*1000/50:3.5*1000/50)),'all');

load('Cnoise_vs_FR0.mat','fit');
Cnoise = (mean(FR_r0)-fit(1))/fit(2);

% back of envelope calculation of upper limit of sum of synaptic
% conductances
gsyn_sum = 0.21;%*(data_FR_masked((best_loc-1)*5 + 1) - steady_state_FR)./STRF_FR;
% gsyn_sum = 0.21*(data_onset - mean(FR_r0))./STRF_onset;


%% Adaptation conductance search

% for 2d adaptation search, use single-channel case
ranges = cell(1,4);
for c = 1:4
    if c == best_loc
        ranges{c} = gsyn_sum;
    else
        ranges{c} = 0; 
    end
end

rmpath(genpath('ICSimStim'));

%% varied parameters

varies = struct;

varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;

varies(end+1).conxn = 'C';
varies(end).param = 'noise';
varies(end).range = Cnoise;

varies(end+1).conxn = 'C';
varies(end).param = 'G_inc';
varies(end).range = 0:0.0004:0.002;
    
varies(end+1).conxn = 'C';
varies(end).param = 'tau_ad';
varies(end).range = 175;

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

adaptFlag = 1;
allFlag = 0;
restrict_vary_flag = 0;

[simdata,DirPart] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,Spks_masked,...
    dataCh,0,folder,subject,'-adaptation-search',best_loc,adaptFlag,allFlag,restrict_vary_flag);

%% compare model and data PSTHs

temp = {simdata.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));

nvaried = size(simdata(mixedIdx).annot,1);

% average PSTH between target identities
mean_PSTH_masked = mean(data_PSTH_masked);
mean_PSTH = mean(data_PSTH);

t_vec = -1:1/50:4;
t_inds = find(t_vec >= onset & t_vec < 3);

startInd = find(t_vec >= 0,1);

% colocated psths

figure('visible','on');
for vv = 1:nvaried
    
    mixedT = 0:0.02:0.02*(length(simdata(mixedIdx).PSTH(vv,:))-1);
    endInd = startInd + size(simdata(mixedIdx).PSTH,2) - 1;
    
    subplot(1,3,2);
    plot(mixedT,scale_data*mean_PSTH_masked(startInd:endInd),'linewidth',1);
    hold on;
    plot(mixedT,scale_model/2*simdata(mixedIdx).PSTH(vv,:),'linewidth',1);
    ylim([0 300]); title('Co-located response');
    legend('Data','Model');
    xlabel('Time (s)'); ylabel('FR'); 
    
    cleanT = 0:0.02:0.02*(length(simdata(targetIdx).PSTH(vv,:))-1);
    endInd = startInd + size(simdata(targetIdx).PSTH,2) - 1;
    
    subplot(1,3,1);
    plot(cleanT,scale_data*mean_PSTH(startInd:endInd),'linewidth',1);
    hold on;
    plot(cleanT,scale_model/2*simdata(targetIdx).PSTH(vv,:),'linewidth',1);
    ylim([0 300]); title('Clean response');
    legend('Data','Model');
    xlabel('Time (s)'); ylabel('FR'); 
    
    annotation('textbox',[.675 .85 .2 .1],...
        'string',simdata(targetIdx).annot(vv,:),...
        'FitBoxToText','on',...
        'LineStyle','none')
    
    saveas(gcf,fullfile(filesep,DirPart,['PSTHs_' num2str(vv) '.tiff']));
    clf;
end
close;

save(fullfile([filesep DirPart],'fitting_results.mat'),'Cnoise',...
    'ICdirPath','ICstruc','gsyn_sum','data_FR','data_FR_masked','data_perf');

%% 4D search
gsyn_range = 0:gsyn_sum/5:gsyn_sum;

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

varies(end+1).conxn = 'C';
varies(end).param = 'G_inc';
varies(end).range = G_inc;

varies(end+1).conxn = 'C';
varies(end).param = 'tau_ad';
varies(end).range = tau_ad;

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
varies(end+1).conxn = 'R->C';
varies(end).param = {'gSYN1','gSYN2','gSYN3','gSYN4'};
varies(end).range = gsyn_sum;

nvaried = {varies(2:end).range};
nvaried = prod(cellfun(@length,nvaried));

allFlag = 0;
restrict_vary_flag = 1;
adaptFlag = 0;
[simdata,DirPart] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,Spks_masked,...
    dataCh,plot_rasters,folder,subject,'-4D-search',[],adaptFlag,allFlag,restrict_vary_flag);

%% performance grids for 4D search

temp = {simdata.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
%maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
perf_clean = []; model_FR = [];

% clean performance and FR
for i = 1:length(targetIdx)
    perf_clean(:,i) = simdata(targetIdx(i)).perf.C;
    model_FR(:,i) = simdata(targetIdx(i)).fr.C;
end

perf_masked = []; model_FR_masked = [];
% mixed performance and FR
if ~isempty(mixedIdx)
for i = 1:length(mixedIdx)
    perf_masked(:,i) = simdata(mixedIdx(i)).perf.C;
    model_FR_masked(:,i) = simdata(mixedIdx(i)).fr.C;
end
end

[~,MSE_clean_perf] = calcModelLoss(perf_clean,data_perf(1:4)');
% calculate loss for FR as a percentage of model FR
[~,MSE_clean_FR] = calcModelLoss(100*(model_FR)./data_FR,100*data_FR./data_FR);

if ~isempty(mixedIdx)
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
    makeParallelPlot(simdata,within_thresh,loss);
else
    best_iterations = zeros(1,sum(within_thresh));
end

for i = 1:length(best_iterations)
    best_iterations(i) = find(loss == temp(i));
end

% make grid of best iterations
makeGrids_bestIteration(simdata,DirPart,data_perf,data_FR,best_iterations,loss);

%% rerun dynasim and obtain rasters for best iteration

% recalculate noise parameter to increase FR and decrease performance
load('Cnoise_vs_FR0.mat','fit');
Cnoise = ((mean(data_FR) - mean(model_FR(best_iterations(1),:)) + mean(FR_r0))-fit(1))/fit(2);

gsyn_strs = cellfun(@str2num,extractAfter({simdata(targetIdx(1)).annot{:,end}},'RC_{gSYN} = '),'UniformOutput',false);
best_gsyns = gsyn_strs{best_iterations(1)};

varies = struct;

varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;

varies(end+1).conxn = 'C';
varies(end).param = 'noise';
varies(end).range = Cnoise;

varies(end+1).conxn = 'C';
varies(end).param = 'G_inc';
varies(end).range = G_inc;

varies(end+1).conxn = 'C';
varies(end).param = 'tau_ad';
varies(end).range = tau_ad;

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

allFlag = 1; % run all spots on grid
restrict_vary_flag = 0;
plot_rasters = 1;
adaptFlag = 0;

[simdata,DirPart] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,...
    Spks_masked,dataCh,plot_rasters,folder,subject,'-best-iteration',[],adaptFlag,allFlag,restrict_vary_flag);

%% make full spatial grid

temp = {simdata.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
%maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
perf_clean = []; model_FR = [];

% clean performance and FR
for i = 1:length(targetIdx)
    perf_clean(:,i) = simdata(targetIdx(i)).perf.C;
    model_FR(:,i) = simdata(targetIdx(i)).fr.C;
end

perf_masked = []; model_FR_masked = [];
% mixed performance and FR
if ~isempty(mixedIdx)
for i = 1:length(mixedIdx)
    perf_masked(:,i) = simdata(mixedIdx(i)).perf.C;
    model_FR_masked(:,i) = simdata(mixedIdx(i)).fr.C;
end
end

[~,MSE_clean_perf] = calcModelLoss(perf_clean,data_perf(1:4)');
% calculate loss for FR as a percentage of model FR
[~,MSE_clean_FR] = calcModelLoss(100*(model_FR)./data_FR,100*data_FR./data_FR);

if ~isempty(mixedIdx)
    [~,MSE_masked_perf] = calcModelLoss(perf_masked,data_perf(5:end)');
    % calculate loss for FR as a percentage of model FR
    [~,MSE_masked_FR] = calcModelLoss(100*(model_FR_masked)./data_FR_masked,100*data_FR_masked./data_FR_masked);
    loss = MSE_clean_perf(:,1) + MSE_clean_FR(:,1) + MSE_masked_perf(:,1) + MSE_masked_FR(:,1);
else
    loss = MSE_clean_perf(:,1) + MSE_clean_FR(:,1);
end

makeGrids(simdata,DirPart,data_perf,data_FR,loss);
