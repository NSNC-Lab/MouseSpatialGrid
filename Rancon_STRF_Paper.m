%0. Load in the stims
LoadStims;

%1. Filter Time Constant eq (Wilmore Paper)
tf = 500 - 105*log10(f_backup);  %Log 10 (see Wilmore Paper)

%2. Calculate the filter length (Rancon Paper)
filter_length = ceil(3*tf(1)) + 1;

%3. Create a matrix with all of the Filters
%Filter Equation (Rancon)
a = exp(-1./tf);
C = 1 - a;

%Onset
filter_bankOn = [];
delta = 1;
w = 0.75;

for k = 1:filter_length
    if k == 1
       filter_bankOn = [filter_bankOn,ones(length(a),3)]; % Basically there is only going to be no place but the start where this function is valid
    else
       filter_bankOn = [filter_bankOn,(-C.*w.*a.^(k-1))'];
    end
end

%Offset
filter_bankOff = [];
delta = 1;
w = 0.75;

for k = 1:filter_length
    if k == 1
       filter_bankOff = [filter_bankOff,-w.*ones(length(a),3)]; % Basically there is only going to be no place but the start where this function is valid
    else
       filter_bankOff = [filter_bankOff,(C.*a.^(k-1))'];
    end
end

%Append delays
filter_bankOff = [filter_bankOff,zeros(length(a),400)]; %40ms delay Taken from Lohse Mgbv STRF
filter_bankOn = [filter_bankOn,zeros(length(a),200)];



%4. Convolve with stimulus represntation

%First index of filter bank is lowest frequency band
%First index of song1 spec should be lowest to highest frequency

convolved_spec_On1 = [];
convolved_spec_Off1 = [];
convolved_spec_On2 = [];
convolved_spec_Off2 = [];

for k = 1:size(song1_spec,2)
    On_Conv = conv(song1_spec(:,k),filter_bankOn(k,:));
    On_Conv = On_Conv(1:end-length(filter_bankOn));

    convolved_spec_On1 = [convolved_spec_On1,On_Conv];%,'same')];
    
    Off_Conv = conv(song1_spec(:,k),filter_bankOff(k,:));
    Off_Conv = Off_Conv(1:end-length(filter_bankOff));

    convolved_spec_Off1 = [convolved_spec_Off1,Off_Conv];%,'same')];

end

for k = 1:size(song2_spec,2)
    On_Conv = conv(song2_spec(:,k),filter_bankOn(k,:));
    On_Conv = On_Conv(1:end-length(filter_bankOn));

    convolved_spec_On2 = [convolved_spec_On2,On_Conv];%,'same')];
    
    Off_Conv = conv(song2_spec(:,k),filter_bankOff(k,:));
    Off_Conv = Off_Conv(1:end-length(filter_bankOff));

    convolved_spec_Off2 = [convolved_spec_Off2,Off_Conv];%,'same')];

end

%5. Sum accross frenecy Bands (Rancon Paper "Linear")



% 
% 
% %1. Sample an exponential filter for excitatory and inhibitory cells
% 
% %Changing filter to be closer to AdapTrans but w/ delays
% 
% filter_length = 10;
% 
% dt = 0.1; %Dynasim Sampling Rate
% dt2 = 5; %Sampling Rate in AdapTrans Paper Cochleograms
% 
% filter_length = filter_length*(dt2/dt);
% 
% filter_l = linspace(0,filter_length-1,filter_length);
% 
% time_length_E = 20;
% time_length_I = 30;
% 
% 
% 
% %t_E_delay = linspace(0,time_length_E,time_length_E/dt);
% %t_I_delay = linspace(0,time_length_I,time_length_I/dt);
% 
% tc = 0.001;
% 
% expo_E = -(exp(filter_l*tc)-1)/(exp((filter_l(end-1))*tc)); %-1 to move down to zero and then scale. This scaling basically acts as an IC gain.
% expo_I = (exp(filter_l*tc)-1)/(exp((filter_l(end-1))*tc));
% 
% expo_E(end) = -sum(expo_E(1:end-1));
% expo_I(end) = -sum(expo_I(1:end-1));
% 
% %Zero pad for delays
% expo_E = [expo_E,zeros(1,time_length_E/dt)];
% expo_I = [expo_I,zeros(1,time_length_I/dt)];
% 
% expo_E = fliplr(expo_E); %The adap trans paper (and many others) do not actually convolve by the definition of convolution so we have to flip to match them.
% expo_I = fliplr(expo_I);
% 
% %Plot the filters
% figure;
% subplot(2,1,1)
% stem(expo_E)
% subplot(2,1,2)
% stem(expo_I)
% 
% %2. convole accross the frequency bands of stimulus
% 
% convolved_spec_On1 = [];
% convolved_spec_Off1 = [];
% convolved_spec_On2 = [];
% convolved_spec_Off2 = [];
% 
% %To propoerly secure casuality we should not use same. The entire filter
% %should slide accross the signal and then on overlap with the signal we
% %should see a sresponse. In order to account for the change in length due
% %to the output of the convolution the other option is to trim the output. I
% %believe as long as the output is zero padded at the signal end this should
% %not be problematic. 
% 
% 
% for k = 1:size(song1_spec,2)
%     On_Conv = conv(song1_spec(:,k),expo_E);
%     On_Conv = On_Conv(1:end-length(expo_E));
% 
%     convolved_spec_On1 = [convolved_spec_On1,On_Conv];%,'same')];
% 
%     Off_Conv = conv(song1_spec(:,k),expo_I);
%     Off_Conv = Off_Conv(1:end-length(expo_I));
% 
%     convolved_spec_Off1 = [convolved_spec_Off1,Off_Conv];%,'same')];
% 
% end
% 
% for k = 1:size(song2_spec,2)
%     On_Conv = conv(song2_spec(:,k),expo_E);
%     On_Conv = On_Conv(1:end-length(expo_E));
% 
%     convolved_spec_On2 = [convolved_spec_On2,On_Conv];%,'same')];
% 
%     Off_Conv = conv(song2_spec(:,k),expo_I);
%     Off_Conv = Off_Conv(1:end-length(expo_I));
% 
%     convolved_spec_Off2 = [convolved_spec_Off2,Off_Conv];%,'same')];
% 
% end


figure;
subplot(2,1,1)
surf(convolved_spec_Off1,'EdgeColor','none')
%x_len = linspace(0,length(song1),length(sum(song1_spec,2)));
%plot(x_len,sum(song1_spec,2));
%xlim([0 length(song1)])
subplot(2,1,2)
surf(convolved_spec_Off2,'EdgeColor','none')
%plot(song1)
%xlim([0 length(song1)])

figure;
subplot(2,1,1)
%surf(convolved_spec_On,'EdgeColor','none')
x_len = linspace(0,length(song1),length(sum(convolved_spec_On1,2)));
plot(x_len,sum(convolved_spec_On1,2));
xlim([0 length(sum(convolved_spec_On1,2))])
subplot(2,1,2)
%surf(convolved_spec_Off,'EdgeColor','none')
plot(song1)
xlim([0 length(song1)])

%3. Get firing rates
convolved_FR_On1 = 80*sum(convolved_spec_On1,2)/max(sum(convolved_spec_On1,2)); %Sum of frequency bands and adjust gain
convolved_FR_Off1 = 80*sum(convolved_spec_Off1,2)/max(sum(convolved_spec_Off1,2));
convolved_FR_On2 = 80*sum(convolved_spec_On2,2)/max(sum(convolved_spec_On2,2)); %Sum of frequency bands and adjust gain
convolved_FR_Off2 = 80*sum(convolved_spec_Off2,2)/max(sum(convolved_spec_Off2,2));

convolved_FR_On1(convolved_FR_On1 < 0) = 0;
convolved_FR_Off1(convolved_FR_Off1<0) = 0;
convolved_FR_On2(convolved_FR_On2 < 0) = 0;
convolved_FR_Off2(convolved_FR_Off2<0) = 0;

figure;
subplot(3,1,1)
plot(convolved_FR_On1)
subplot(3,1,2)
plot(convolved_FR_Off1)
subplot(3,1,3)
plot(song1)

fr_target_on{1} = convolved_FR_On1;
fr_target_off{1} = convolved_FR_Off1;
fr_target_on{2} = convolved_FR_On2;
fr_target_off{2} = convolved_FR_Off2;


%%%%%



save('Rancon_STRF_Paper.mat','fr_target_on','fr_target_off','paramH','paramG','strfGain');

toc