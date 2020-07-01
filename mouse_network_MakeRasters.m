clear
close all

addpath('mechs')
addpath('dependencies')
addpath('eval_scripts')
addpath('genlib')
addpath(genpath('dynasim'))

simFolder = uigetdir('run');

files = dir(simFolder); files(~[files.isdir] | contains({files.name},'.')) = [];
addpath(genpath([simFolder filesep]));
load([simFolder filesep files(1).name filesep 'studyinfo.mat']);

numSims = length(studyinfo.simulations);
varies = studyinfo.base_simulator_options.vary;
clearvars studyinfo

set(0, 'DefaultFigureVisible', 'off')
h = figure('Position',[50,50,850,690]);

for i = 1:length(files)   % for each configuration
    mouse_network_rasteronly([simFolder filesep files(i).name],numSims,varies)
end