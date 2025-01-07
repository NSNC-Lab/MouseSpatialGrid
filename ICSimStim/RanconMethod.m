%%% Inputs

close all

%Decleare gammatone fileter bank
fs = 100000;
fs2 = 64000;
num_filts = 128;
gammaFiltBank = gammatoneFilterBank([0,fs2/2],num_filts,fs2);

%Declare toy stimulus (From linden paper)
Start_pause = zeros(1,(200/1000)*fs);
stimuli = wgn(1,(150/1000)*fs,0);
stimuli_2 = wgn(1,(5/1000)*fs,0);
linear_kernal = linspace(0,1,length(stimuli_2));

stimuli_on = linear_kernal.*stimuli_2;
stimuli_off = fliplr(linear_kernal).*stimuli_2;
silence = zeros(1,(150/1000)*fs);
full_stimuli = [Start_pause,stimuli_on,stimuli,stimuli_off,silence];

figure;
plot(full_stimuli)

%Apply the gammatone filter bank to the stim
audioOut = squeeze(gammaFiltBank(full_stimuli));

%Calculate time and frequency vectors
samplesPerFrame = round(0.005*fs);
samplesOverlap = round(0.0025*fs);

buff = dsp.AsyncBuffer(numel(full_stimuli));
write(buff,squeeze(audioOut)'.^2);

sink = dsp.AsyncBuffer(numel(full_stimuli));

while buff.NumUnreadSamples > 0
   currentFrame = read(buff,samplesPerFrame,samplesOverlap);
   write(sink,mean(currentFrame,1));
end

gammatoneSpec = read(sink);
D = 20*log10(gammatoneSpec');
D = gammatoneSpec';

timeVector = ((samplesPerFrame-samplesOverlap)/fs)*(0:size(D,2)-1);
cf = getCenterFrequencies(gammaFiltBank)./1e3;


% 
% timeVector = linspace(0,length(full_stimuli)/fs,size(audioOut,2));
% cf = getCenterFrequencies(gammaFiltBank)./1e3;


%Display Cochleogram
figure;
surf(timeVector,cf,D,EdgeColor="none")
view([0 90])
xlabel("Time (s)")
ylabel("Frequency (kHz)")





%%% Adapt Trans Filters

%Filter size = 10
%w = 0.5
%a = 0.6
%Hon = delt - Con*w*sum_term_on;
%Hoff = -w*delt+Coff*sum_term_off;

Hon = fliplr([-0.001,-0.005,-0.01,-0.02,-0.04,-0.07,-0.11,-0.18,-0.24,1]);
Hoff = fliplr([0.001,0.005,0.01,0.04,0.07,0.11,0.18,0.24,0.4,-0.5]);






%%% Convolve with freuqncy bands

D_conv_on = [];
D_conv_off = [];

for k = 1:num_filts
    D_conv_on = [D_conv_on;conv(D_conv_strf(k,:),Hon)];
    D_conv_off = [D_conv_off;conv(D_conv_strf(k,:),Hoff)];
end

timeVector2 = ((samplesPerFrame-samplesOverlap)/fs)*(0:size(D_conv_on,2)-1);

%%% Display onset and offset spectrograms

%Onset
figure;
surf(timeVector2,cf,D_conv_on,EdgeColor="none")
view([0 90])
title('Onset')
xlabel("Time (s)")
ylabel("Frequency (kHz)")

%Offset
figure;
surf(timeVector2,cf,D_conv_off,EdgeColor="none")
view([0 90])
title('Offset')
xlabel("Time (s)")
ylabel("Frequency (kHz)")



%%% Sum accross frequency bands and show 2D response
D_conv_off_2D = sum(D_conv_off);
D_conv_on_2D = sum(D_conv_on);


figure;
plot(timeVector2,D_conv_on_2D)
xlabel("Time (s)")
ylabel("Sum across frequencies")
title('Summation plot ON')

figure;
plot(timeVector2,D_conv_off_2D)
xlabel("Time (s)")
ylabel("Sum across frequencies")
title('Summation plot OFF')



%%% Convolve STRF with frequency bands

tic

if ~contains(pwd,'ICSimStim'), cd('ICSimStim'); end
%clearvars
%clc
close all

addpath(genpath('strflab_v1.45'))
addpath('../genlib')
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

paramH.alpha = 0.0055; % time constant of temporal kernel [s] 0.0097
paramH.N1 = 5;
paramH.N2 = 7;
paramH.SC1 = 1; %Nominally 1
paramH.SC2 = 1.2;  % increase -> more inhibition %0.88 in paper
strfGain = 0.1;

% frequency parameters - static
paramG.BW = 2000;  % Hz
paramG.BSM = 5.00E-05; % 1/Hz=s best spectral modulation
paramG.f0 = 4300; % ~strf.f(30)

% make STRF
strf = STRFgen_V2(paramH,paramG,cf,0.00015);
figure;
plot(sum(strf.w1))


D_on_conv_strf = [];
D_off_conv_strf = [];

for k = 1:num_filts
    D_on_conv_strf = [D_on_conv_strf;conv(D_conv_on(k,:),sum(strf.w1))];
    D_off_conv_strf = [D_off_conv_strf;conv(D_conv_off(k,:),sum(strf.w1))];
end

timeVector2 = ((samplesPerFrame-samplesOverlap)/fs)*(0:size(D_off_conv_strf,2)-1);

figure;
surf(timeVector2,cf,D_on_conv_strf,EdgeColor="none")
view([0 90])
title('Onset')
xlabel("Time (s)")
ylabel("Frequency (kHz)")

figure;
surf(timeVector2,cf,D_off_conv_strf,EdgeColor="none")
view([0 90])
title('Offset')
xlabel("Time (s)")
ylabel("Frequency (kHz)")


