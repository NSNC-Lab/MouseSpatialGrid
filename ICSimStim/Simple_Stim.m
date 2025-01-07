tic

if ~contains(pwd,'ICSimStim'), cd('ICSimStim'); end

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

paramH.alpha = 0.0042; % time constant of temporal kernel [s] 0.0097
paramH.N1 = 5;
paramH.N2 = 7;
paramH.SC1 = 1; %Nominally 1
paramH.SC2 = 1;  % increase -> more inhibition %0.88 in paper

strfGain = 0.1/50;

% frequency parameters - static
paramG.BW = 32000;  % Hz
paramG.BSM = 8.0E-06; % 1/Hz=s best spectral modulation
paramG.f0 = 32000; % ~strf.f(30)

%IB 12/2
%Updated for broadband STRF according to Olsen & HausenStaub
%paramG.BW = 16000;  % Hz

%They say their sampling rate is 192000

fs = 192000;

%200 ms starting pause
Start_pause = zeros(1,(200/1000)*fs);


stimuli = wgn(1,(150/1000)*fs,0);

%Ramp on
stimuli_2 = wgn(1,(5/1000)*fs,0);
linear_kernal = linspace(0,1,length(stimuli_2));

stimuli_on = linear_kernal.*stimuli_2;
stimuli_off = fliplr(linear_kernal).*stimuli_2;



%Ramp off

silence = zeros(1,(150/1000)*fs);

full_stimuli = [Start_pause,stimuli_on,stimuli,stimuli_off,silence,stimuli_on,stimuli,stimuli_off,silence,stimuli_on,stimuli,stimuli_off,silence,stimuli_on,stimuli,stimuli_off,silence,Start_pause];

full_stimuli = [full_stimuli,full_stimuli,full_stimuli];

[full_stimuli_spec,t,f]=STRFspectrogram_Olsen(full_stimuli/rms(full_stimuli)*targetlvl,fs);

%[full_stimuli_spec,t,f] = spectrogram(full_stimuli,32,16);
