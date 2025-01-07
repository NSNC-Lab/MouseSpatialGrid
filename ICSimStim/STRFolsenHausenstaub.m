tic

if ~contains(pwd,'ICSimStim'), cd('ICSimStim'); end
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

[full_stimuli_spec,t,f]=STRFspectrogram_Olsen(full_stimuli/rms(full_stimuli)*targetlvl,fs);


figure;
surf(f,t,full_stimuli_spec,'EdgeColor','none')


specs.songs{1} = full_stimuli_spec;
specs.dims = size(full_stimuli_spec);
specs.t = t;
specs.f = f;
save('preprocessed_stim_Olsen_200k.mat','specs')
% end




% make STRF
strf = STRFgen_V2(paramH,paramG,specs.f,specs.t(2)-specs.t(1));
strf2 = STRFgen_V2_flipped(paramH,paramG,specs.f,specs.t(2)-specs.t(1));

strf.w1 = strf.w1*strfGain;
strf2.w1 = strf2.w1*strfGain;

figure('Position',[100,100,500,500]);
plot(1000*strf.t,sum(strf.w1)); hold on;
xlabel('Time (ms)'); ylabel('Weight'); title('MGB STRF (Figure 3c)')
savefig(gcf,'MGB STRF Figure 3c, v3');

figure('Position',[100,100,500,500]);
plot(1000*strf2.t,sum(strf2.w1)); hold on;
xlabel('Time (ms)'); ylabel('Weight'); title('MGB STRF (Figure 3c)')
savefig(gcf,'MGB STRF Figure 3c, v3');

paramSpk.t_ref = 1.5;
paramSpk.t_ref_rel = 0.5;
paramSpk.rec = 4;

%% Run simulation script
mean_rate = 0.1;

tuningParam.strf = strf;
tuningParam.type = tuning;
tuningParam.sigma = sigma;

%[~,~,fr_target_on{1},fr_target_off{1}] = STRFconvolve_V2(strf,specs.songs{1}*stimGain,mean_rate,1,[],paramSpk.t_ref,paramSpk.t_ref_rel,paramSpk.rec);
%[~,~,fr_target_on{2},fr_target_off{2}] = STRFconvolve_V2(strf,specs.songs{2}*stimGain,mean_rate,1,[],paramSpk.t_ref,paramSpk.t_ref_rel,paramSpk.rec);
[fr_target_on{1},fr_target_off{1}] = STRFconvolve_V2(strf,specs.songs{1}*stimGain,mean_rate);

[fr_target_on2{1},fr_target_off2{1}] = STRFconvolve_V2(strf2,specs.songs{1}*stimGain,mean_rate);



figure('position',[100,100,500,500])
subplot(5,1,1)
plot(full_stimuli)
subplot(5,1,2)
plot(fr_target_on{1})
subplot(5,1,3)
plot(fr_target_off{1})
subplot(5,1,4)
plot(fr_target_on2{1})
subplot(5,1,5)
plot(fr_target_off2{1})


figure('position',[100,100,500,500])

full_0_idx = [];

for mm = 1:length(full_stimuli)-1
    if (full_stimuli(mm+1) ~= 0 && full_stimuli(mm) == 0) || (full_stimuli(mm+1) == 0 && full_stimuli(mm) ~= 0)
        full_0_idx = [full_0_idx,mm];
    end
end


on_0_idx = [];

for mm = 1:length(fr_target_on{1})-1
    if (fr_target_on{1}(mm+1) > 0 && fr_target_on{1}(mm) == 0) || (fr_target_on{1}(mm+1) == 0 && fr_target_on{1}(mm) > 0)
        on_0_idx = [on_0_idx,mm];
    end
end

on2_0_idx = [];

for mm = 1:length(fr_target_on2{1})-1
    if (fr_target_on2{1}(mm+1) > 0 && fr_target_on2{1}(mm) == 0) || (fr_target_on2{1}(mm+1) == 0 && fr_target_on2{1}(mm) > 0)
        on2_0_idx = [on2_0_idx,mm];
    end
end


subplot(3,1,1)
plot(full_stimuli); hold on

for kk = 1:length(full_0_idx)
    plot([full_0_idx(kk),full_0_idx(kk)],[-5,5],'k--'); hold on
end

xlim([1 314880])
subplot(3,1,2)
plot(fr_target_on{1}); hold on

for kk = 1:length(on_0_idx)/2
    plot([on_0_idx(kk*2-1),on_0_idx(kk*2-1)],[0,20],'k--'); hold on
end

xlim([1 16573])
subplot(3,1,3)
plot(fr_target_on2{1}); hold on

for kk = 1:length(on2_0_idx)/2
    plot([on2_0_idx(kk*2-1),on2_0_idx(kk*2-1)],[0,20],'k--'); hold on
end

xlim([1 16573])


save('default_Olsen_STRF_with_offset_200k.mat','fr_target_on','fr_target_off','paramH','paramG','strfGain');

toc