%Add a path to the GA functions
addpath('GA_Functions\')
clc;
clear all

%Create a function to handle all of this stuff

%Move to the correct simulation path
cd('..');
warning('off','all');

% change current directory to folder where this script is stored
mfileinfo = mfilename('fullpath');
mfiledir = strsplit(mfileinfo,filesep);
% cd(fullfile(mfiledir{1:end-1}));

dynasimPath = 'DynaSim-master';

addpath('mechs');
addpath('resampled-stimuli');
addpath(genpath('ICSimStim'));
addpath('genlib');
addpath(genpath(dynasimPath));
addpath('cSPIKE'); InitializecSPIKE;
addpath('plotting');
addpath('subfunctions');
addpath('MICalculation');

plot_all = 0;
%Uses "masked" tuning curves instead of all clean tuning curves as input
toggle_real = 0;