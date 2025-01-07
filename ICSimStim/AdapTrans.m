
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

%full_stimuli = [Start_pause,stimuli_on,stimuli,stimuli_off,silence,stimuli_on,stimuli,stimuli_off,silence,stimuli_on,stimuli,stimuli_off,silence,stimuli_on,stimuli,stimuli_off,silence,Start_pause];

%full_stimuli = bandpass(full_stimuli,[4000,64000],fs);

tau1 = 0.006;
tau2 = 0.01;

x = linspace(-0.1,0,1000);
x2 = linspace(-0.1,0,1000);

k1 = exp((x)/tau1);
k2 = exp((x2)/tau2)*0.5;

k1 = [k1];
k2 = [k2];

figure;
subplot(2,1,1)
plot(k1)
subplot(2,1,2)
plot(k2)

%Zero padd k1 to be the same length a k2
%dif_len = length(k2)-length(k1);

%k1 = [k1,zeros(1,dif_len)];


%figure;
%plot(k1)

%figure;
%plot(k2)



%\/\/\/\/\/\/\/\/\/\/
% Linden Method
%\/\/\/\/\/\/\/\/\/\/

%1. Take the sound envelope s(t)
%env_full_stimuli = envelope(full_stimuli,1000,'rms');

%Try toy stim for now
env_full_stimuli = [ones(1,10000)];

%2. Intensity Integration
rI = conv(env_full_stimuli,k1);

%3. Adaptive Gain
%gt = 1./(1+conv(rI,k2));
gt2 = conv(rI,k1);

gt2_temp  = 1./(max(rI)+gt2(1000:11998));

%4. Get Intensity Gain control
rIA = gt2_temp.*rI;

figure
subplot(6,1,1)
plot(k1)
subplot(6,1,2)
plot(k2)
subplot(6,1,3)
plot(env_full_stimuli)
subplot(6,1,4)
plot(rI)
subplot(6,1,5)
plot(gt2_temp)
%xlim([1000,12000])
subplot(6,1,6)
plot(rIA)




%Notes: Window length must be the same. Trying 'same' comvolution to keep
%window size