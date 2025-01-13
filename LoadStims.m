% Simulate spike times based on a model neuron with a (1) STRF and a (2)
% spatial tuning curve. I.e. this neuron has both spatial and
% spectral-temporal tuning.

% note:
% large bottleneck lies in r/w to network drive
tic

if ~contains(pwd,'ICSimStim'), cd('ICSimStim'); end
%clearvars
%clc
close all

addpath(genpath('strflab_v1.45'))
addpath('../genlib')
% addpath('../fixed-stimuli')
addpath('../resampled-stimuli')
addpath('../plotting')
addpath('../cSPIKE')
InitializecSPIKE;
dataSaveLoc = pwd; %local save location

% Spatial tuning curve parameters
sigma = 24; %60 for bird but 38 for mouse
tuning = 'mouse'; %'bird' or 'mouse'
stimGain = 0.5;
targetlvl = 0.01;
maskerlvl = 0.01; %default is 0.01
maxWeight = 1; %maximum mixed tuning weight; capped at this level.

paramH.alpha = 0.0097; % time constant of temporal kernel [s] 0.0097
paramH.N1 = 5;
paramH.N2 = 7;
paramH.SC1 = 1; %Nominally 1
paramH.SC2 = 0.88;  % increase -> more inhibition %0.88 in paper

strfGain = 0.1;

% frequency parameters - static
paramG.BW = 2000;  % Hz
paramG.BSM = 5.00E-05; % 1/Hz=s best spectral modulation
paramG.f0 = 4300; % ~strf.f(30)
    
[song1,~] = audioread('200k_target1.wav');
[song2,~] = audioread('200k_target2.wav');

for trial = 1:10
    [masker,fs] = audioread(['200k_masker' num2str(trial) '.wav']);
    [spec,~,~] = STRFspectrogram(masker/rms(masker)*maskerlvl,fs);
    masker_specs{trial} = spec;
end

[song1_spec,t_backup,f_backup]=STRFspectrogram(song1/rms(song1)*targetlvl,fs);
[song2_spec,~,~]=STRFspectrogram(song2/rms(song2)*targetlvl,fs);
specs.songs{1} = song1_spec;
specs.songs{2} = song2_spec;
specs.maskers = masker_specs;
specs.dims = size(song1_spec);
specs.t = t_backup;
specs.f = f_backup;
save('preprocessed_stim_200k.mat','specs')
