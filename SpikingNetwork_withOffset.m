%% Initialize

mfileinfo = mfilename('fullpath');
mfiledir = fileparts(mfileinfo);
% cd(mfiledir);

dynasimPath = '../DynaSim';

addpath('mechs');
addpath('resampled-stimuli');
addpath(genpath('ICSimStim'));
addpath('genlib');
addpath(genpath(dynasimPath));
addpath('cSPIKE'); InitializecSPIKE;
addpath('plotting');
addpath('subfunctions');

load('default_STRF_with_offset_200k.mat');
dt = 0.2; %ms

% what if we wanted to downsample the firing rates?
if dt ~= 0.1
    dsamp_fac = dt/0.1;
    for m = 1:10
        fr_masker{m} = downsample(fr_masker{m},dsamp_fac);
    end
    for t = 1:2
    fr_target_on{t} = downsample(fr_target_on{t},dsamp_fac);
    fr_target_off{t} = downsample(fr_target_off{t},dsamp_fac);
    end
end

% study_dir: folder under 'run' where m files and input spikes for simulations are written and saved
study_dir = fullfile(pwd,'run','3-channel');

if exist(study_dir, 'dir'), msg = rmdir(study_dir, 's'); end
mkdir(fullfile(study_dir, 'solve'));

% expName: folder under 'simData' where results are saved
expName = '12-01-23 3-channel test V2';

newStrfGain = strfGain;
simDataDir = [pwd filesep 'simData' filesep expName];

if ~exist(simDataDir,'dir'), mkdir(simDataDir); end

%% Run .m file to generate options and varies structs for simulations
addpath('params-3-channel');
params_3channel_Full;

% for figures in paper

% addpath('params');
%Figure 4a, params_4a;
%Figure 4b, params_4b;
% Figure 4c, params_4c;

% For Figure 5, params_5
% For Figure 6, params_6
% For Figure 7, params_7
% for Figure 8, params_8

% for masked params i don't know what to do with
% params_MaskedPerf;
% params_Masked_varyOnsetPV;
% params_Masked_varyOnsetNoise;
% params_DepressiveStr;
% params_ExpFig3;

%% create input spikes from STRFs
% concatenate spike-time matrices, save to study dir

prepInputData;

%% run simulation

options.strfGain = newStrfGain;
options.dt = dt;

if isempty(options.locNum), options.time_end = size(spks,1)*dt; %ms;
else, options.time_end = padToTime * numel(options.locNum); end
[snn_out,s] = columnNetwork_V2(study_dir,varies,options,netcons);

%% post-process for performance and firing results

postProcessSims;
