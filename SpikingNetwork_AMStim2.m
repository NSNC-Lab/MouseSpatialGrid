% This version will help you create the AM stimulus responses in Figure 6
% in the bioRxiv paper. - Jio

%This (2) script is just used for testing

%% Initialize

% change current directory to folder where this script is stored
mfileinfo = mfilename('fullpath');
mfiledir = strsplit(mfileinfo,filesep);
% cd(fullfile(mfiledir{1:end-1}));

dynasimPath = ['..' filesep 'DynaSim'];

addpath('mechs'); addpath(genpath('ICSimStim'));
addpath('genlib'); addpath('plotting'); addpath(genpath(dynasimPath));
addpath('cSPIKE'); InitializecSPIKE;
addpath('plotting');
addpath('params-AM');
addpath('AM-code');
addpath('subfunctions');

options.dt = 0.1; %ms

% study_dir: folder under 'run' where m files and input spikes for simulations are written and saved
study_dir = fullfile(pwd,'run','single-channel-AM-stim');

if exist(study_dir, 'dir'), msg = rmdir(study_dir, 's'); end
mkdir(fullfile(study_dir, 'solve'));

load('AM_stim_FR_traces.mat');
load('peakDrv_AM.mat');
load('AM_sigs.mat')

% expName: folder under 'simData' where results are saved
expName = '08-30-23 fix dynamics, vary strength';

% for newStrfGain = strfGains
% simDataDir = [pwd filesep 'simData' filesep expName ' ' num2str(newStrfGain)];

strfGain = 0.08; %strfGain;
simDataDir = [pwd filesep 'simData-AM' filesep expName];

if ~exist(simDataDir,'dir'), mkdir(simDataDir); end

padToTime = (t_stim + 0.250 + 0.250) * 1000; %ms, first 250ms accts for padding in STRF_AMStim, second 250ms is for padding at the end

if length(fr_target_on) > 2 % AM stim with steady modulation
    nStim = 5;
    nTrials = 100;
    TrialsPerStim = nTrials/5;
else
    nStim = 2; % irregular AM stim
    nTrials = 20;
    TrialsPerStim = nTrials/2;
end

%% Run .m file to generate options and varies structs for simulations

% params_AM_adjustPVDynamics;
% params_AM_noDep;
params_AM_best_onoff_OnOff;
% params_AM_varyingStrengths;

%% create input spikes from STRFs
% concatenate spike-time matrices, save to study dir

padToTime = 3500; % [ms]

% ICfiles.mat contains names of spatial grid configs: s[targetloc]m[maskerloc]
% See 'config_idx_reference.JPG' for indexes

% options.locNum is defined in params .m file
% if it's empty, default to running all 24 configs (including masker-only
% trials)

if ~isempty(options.locNum)
    subz = options.locNum;
else
    subz = 1:24;
end

% concatenate spike-time matrices, save to study_dir
prepInputData_AMStim;

%% run simulation

options.strfGain = strfGain;
options.dt = dt;
options.TrialsPerStim = TrialsPerStim;

if isempty(options.locNum), options.time_end = size(spks,1)*dt; %ms;
else, options.time_end = padToTime; end
[snn_out,s] = columnNetwork_paper_onoff(study_dir,varies,options,netcons);

%% post-process for performance and firing results

options.strfGain = strfGain; % store in options for logging
postProcessAM;
