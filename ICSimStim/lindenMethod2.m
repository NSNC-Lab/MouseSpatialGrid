tau1 = 0.006;
tau2 = 0.01;

x = linspace(-0.1,0,1000);
x2 = linspace(-0.1,0,1000);

k1 = exp((x)/tau1);
k2 = exp((x2)/tau2)*0.5;

k1 = [k1];
k2 = [k2];

%\/\/\/\/\/\/\/\/\/\/
% Linden Method
%\/\/\/\/\/\/\/\/\/\/

%1. Take the sound envelope s(t)
%env_full_stimuli = envelope(full_stimuli,1000,'rms');

%Try toy stim for now
env_full_stimuli = ones(1,10000);

%2. Intensity Integration
rI = conv(env_full_stimuli,k1);
rI = (1/max(rI)).*rI;

%3. Adaptive Gain
%gt = 1./(1+conv(rI,k2));
gt2_temp = conv(rI,k1);
max_gt2 = max(gt2_temp);
gt2_temp = gt2_temp(length(x2):length(env_full_stimuli)+length(x2)+length(x)-2);
gt2_temp = (-gt2_temp)+max_gt2;

%gt2_temp  = 1./(1+gt2(length(x2):length(env_full_stimuli)+length(x2)+length(x)-2));

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