%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>%

%Create Adapt Trans filters that are on the same timescale as the STRFS.

%STRFS are currently on the sampling rate of 0.1 ms. 3.5 second timescale
%translated to a STRF that is 35000 samples long.
%Time to peak for excitatory cells is around 20ms 
%Time to peak for inhibitory cells is ______

%For now we are not worried about the frequency componnet so we will do
%frequency band wise convolution with the stimuli and then sum accross the
%frequency bands

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>%

%1. Sample an exponential filter for excitatory and inhibitory cells

time_length_E = 20;
time_length_I = 30;

dt = 0.1; %Sampling Rate

t_E = linspace(0,time_length_E,time_length_E/dt);
t_I = linspace(0,time_length_I,time_length_I/dt);

expo_E = -(exp(t_E*0.1)-1)/(exp((time_length_E-dt)*0.1)); %-1 to move down to zero and then scale. This scaling basically acts as an IC gain.
expo_I = (exp(t_I*0.1)-1)/(exp((time_length_I-dt)*0.1));

expo_E(end) = -sum(expo_E(1:end-1));
expo_I(end) = -sum(expo_I(1:end-1));

expo_E = fliplr(expo_E); %The adap trans paper (and many others) do not actually convolve by the definition of convolution so we have to flip to match them.
expo_I = fliplr(expo_I);

%Plot the filters
figure;
subplot(2,1,1)
stem(t_E,expo_E)
subplot(2,1,2)
stem(t_I,expo_I)

%2. convole accross the frequency bands of stimulus
%Simple_Stim; %Get the simple stimulus
Paper_Stim;

convolved_spec_On = [];
convolved_spec_Off = [];

%To propoerly secure casuality we should not use same. The entire filter
%should slide accross the signal and then on overlap with the signal we
%should see a sresponse. In order to account for the change in length due
%to the output of the convolution the other option is to trim the output. I
%believe as long as the output is zero padded at the signal end this should
%not be problematic. 


for k = 1:size(full_stimuli_spec,2)
    On_Conv = conv(full_stimuli_spec(:,k),expo_E);
    On_Conv = On_Conv(1:end-length(expo_E));

    convolved_spec_On = [convolved_spec_On,On_Conv];%,'same')];
    
    Off_Conv = conv(full_stimuli_spec(:,k),expo_I);
    Off_Conv = Off_Conv(1:end-length(expo_I));

    convolved_spec_Off = [convolved_spec_Off,Off_Conv];%,'same')];

end


figure;
subplot(2,1,1)
%surf(convolved_spec_On,'EdgeColor','none')
x_len = linspace(0,length(full_stimuli),length(sum(full_stimuli_spec,2)));
plot(x_len,sum(full_stimuli_spec,2));
xlim([0 length(full_stimuli)])
subplot(2,1,2)
%surf(convolved_spec_Off,'EdgeColor','none')
plot(full_stimuli)
xlim([0 length(full_stimuli)])

figure;
subplot(2,1,1)
%surf(convolved_spec_On,'EdgeColor','none')
x_len = linspace(0,length(full_stimuli),length(sum(convolved_spec_On,2)));
plot(x_len,sum(convolved_spec_On,2));
xlim([0 length(sum(convolved_spec_On,2))])
subplot(2,1,2)
%surf(convolved_spec_Off,'EdgeColor','none')
plot(full_stimuli)
xlim([0 length(full_stimuli)])

%3. Get firing rates
convolved_FR_On = 300*sum(convolved_spec_On,2)/max(sum(convolved_spec_On,2)); %Sum of frequency bands and adjust gain
convolved_FR_Off = 300*sum(convolved_spec_Off,2)/max(sum(convolved_spec_Off,2));

convolved_FR_On(convolved_FR_On < 0) = 0;
convolved_FR_Off(convolved_FR_Off<0) = 0;

figure;
subplot(2,1,1)
plot(convolved_FR_On)
subplot(2,1,2)
plot(convolved_FR_Off)

fr_target_on = {convolved_FR_On};
fr_target_off = {convolved_FR_Off};

save('SimpleStimInput.mat','fr_target_on','fr_target_off','paramH','paramG','strfGain');