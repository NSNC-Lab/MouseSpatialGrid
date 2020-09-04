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
folder = ['Data-fitting' filesep subject];

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
data_FR_masked = zeros(1,16);
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
load([ICdirPath ICstruc(ICs(best_loc)).name],'avgSpkRate')
STRF_FR_clean = avgSpkRate;

% collocated STRF rates
ICs = find((cellfun(@(x) x(2),{ICstruc.name}) == cellfun(@(x) x(4),{ICstruc.name})));
load([ICdirPath ICstruc(ICs(best_loc)).name],'avgSpkRate')
STRF_FR_coloc = avgSpkRate;

frac_clean = (data_FR(best_loc) - mean(FR_r0)) / mean(STRF_FR_clean);
frac_coloc = (data_FR_masked(1+5*(best_loc-1)) - mean(FR_r0)) / mean(STRF_FR_coloc);

gsyn_sum = max([frac_coloc*0.18 frac_clean*0.18]);

%% Estimate cortical noise for co-located trials

varies = struct;

ranges = cell(1,4);
for c = 1:4
    if c == best_loc
        ranges{c} = 0.18*(data_FR(best_loc)-mean(FR_r0))/mean(STRF_FR_clean);
    else
        ranges{c} = 0;
    end
end

varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;

varies(end+1).conxn = 'C';
varies(end).param = 'noise';
varies(end).range = Cnoise + [0:0.2:2];

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

subz = find(contains({ICstruc.name},['s' num2str(best_loc) 'm' num2str(best_loc)])); % co-located cases

[simdata] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,...
    Spks_masked,dataCh,1,folder,'-collocated-noise',subz,[],0);

% Find best colocated cortical noise

temp = {simdata.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content
colocIdx = find((cellfun(@(x) x(2),temp) == cellfun(@(x) x(4),temp)));

[~,ind] = min(abs(simdata(colocIdx).fr.C - data_FR_masked(1 + 5*(best_loc-1))));
Cnoise2 = varies(2).range(ind);

% save parameters and data for manual fitting

save(fullfile(folder,'default_parameters.mat'),'gain','ICdir','Cnoise','Cnoise2','gsyn_sum',...
    'data_FR','data_FR_masked','data_perf','best_loc','dataCh');
