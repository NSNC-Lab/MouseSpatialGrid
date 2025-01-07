load('all_units_info_updated_manual.mat','all_data')

%77 = 327L CH 18
%78 = 327L CH 20
a = all_data(78).ctrl_tar1_masked_timestamps;
b = all_data(78).ctrl_tar1_timestamps;
nbins = 125;

%  this counts the number of spikes in 3 second, you need to devide it by 3
%  to get the FR
[count_ctrl_tar2_clean, ~, ~] = calc_num_spikes(a(:,:,4),nbins);
% [count_ctrl_tar2_clean, ~, ~] = calc_num_spikes(b,nbins);
FR_ctrl_tar2_clean = count_ctrl_tar2_clean/3;

degree_list_x = [90 45 0 -90];
numBootstraps = 500;

%Load in comparisons
state1 = load('327_CH20_AbstractRun_Noramalized_Fit_Lam4_20_Lam5_30.mat');
state2 = load('327_CH20_AbstractRun_Noramalized_Fit.mat');
state3 = load('327_CH20_AbstractRun_Noramalized_Fit_Lam5_5.mat');
state4 = load('327_CH20_AbstractRun_Noramalized_Fit_Lam4_5_Lam5_20.mat');


% FRS1 = state1.state.frGridHistory{end}(1,:);
% FRS2 = state2.state.frGridHistory{end}(1,:);
% FRS3 = state3.state.frGridHistory{end}(1,:);
% FRS4 = state4.state.frGridHistory{end}(1,:);

FRS1 = state1.state.frGridHistory{end}(2,:);
FRS2 = state2.state.frGridHistory{end}(2,:);
FRS3 = state3.state.frGridHistory{end}(2,:);
FRS4 = state4.state.frGridHistory{end}(2,:);


%Find grid w/ best fitness


% you need to modify the code "Calc_cen_mod_ERRF_bootstrap_raw"  by adding  the variable called "bootstrapMeans" as its output like this:
% function [errf_width_ctrl_all, mod_depth_ctrl_all, centroid_ctrl_all, bootstrapMeans]= calc_cen_mod_ERRF_bootstrap_raw(numBootstraps, counts,degree_list_x)  
[~, ~, ~, samp_mean]= calc_cen_mod_ERRF_bootstrap_raw(numBootstraps, FR_ctrl_tar2_clean ,degree_list_x);  

[~, ~, ~, samp_mean1]= calc_cen_mod_ERRF_bootstrap_raw(numBootstraps, FRS1 ,degree_list_x);  

[~, ~, ~, samp_mean2]= calc_cen_mod_ERRF_bootstrap_raw(numBootstraps, FRS2 ,degree_list_x);  

[~, ~, ~, samp_mean3]= calc_cen_mod_ERRF_bootstrap_raw(numBootstraps, FRS3 ,degree_list_x);  

[~, ~, ~, samp_mean4]= calc_cen_mod_ERRF_bootstrap_raw(numBootstraps, FRS4 ,degree_list_x);  


normalizer_real_data = max(mean(samp_mean));
normalizer_real_data1 = max(mean(samp_mean1));
normalizer_real_data2 = max(mean(samp_mean2));
normalizer_real_data3 = max(mean(samp_mean3));
normalizer_real_data4 = max(mean(samp_mean4));

samp_mean = samp_mean./normalizer_real_data;
samp_mean1 = samp_mean1./normalizer_real_data1;
samp_mean2 = samp_mean2./normalizer_real_data2;
samp_mean3 = samp_mean3./normalizer_real_data3;
samp_mean4 = samp_mean4./normalizer_real_data4;

% then you can then use "Calc_CI" to get the confidence interval at each speaker location.
std_clean_speaker1 = std(samp_mean(:,1));
std_clean_speaker2 = std(samp_mean(:,2));
std_clean_speaker3 = std(samp_mean(:,3));
std_clean_speaker4 = std(samp_mean(:,4));
%CI_clean_speaker4 = calc_CI(mean(:,4), 0.95);

% For masked condition, you need to know that "a = all_data(57).ctrl_tar1_masked_timestamps;"
% has a 10*4*4 structure, so the first 4 is target location and the second
% 4 is the masker location, they all follow the same order, 90 45 0 and -90

% so you could do this:
[count_ctrl_tar1_makser90, ~, ~] = calc_num_spikes(a(:,:,1),nbins);
FR_ctrl_tar1_masker90 = count_ctrl_tar1_makser90/3;

% Example data
x = 1:4;                   % X-axis values
y = [mean(samp_mean(:,1)), mean(samp_mean(:,2)), mean(samp_mean(:,3)), mean(samp_mean(:,4))];  % Y-axis values
y1 = [mean(samp_mean1(:,1)), mean(samp_mean1(:,2)), mean(samp_mean1(:,3)), mean(samp_mean1(:,4))];  % Y-axis values
y2 = [mean(samp_mean2(:,1)), mean(samp_mean2(:,2)), mean(samp_mean2(:,3)), mean(samp_mean2(:,4))];  % Y-axis values
y3 = [mean(samp_mean3(:,1)), mean(samp_mean3(:,2)), mean(samp_mean3(:,3)), mean(samp_mean3(:,4))];  % Y-axis values
y4 = [mean(samp_mean4(:,1)), mean(samp_mean4(:,2)), mean(samp_mean4(:,3)), mean(samp_mean4(:,4))];  % Y-axis values

errors = [std_clean_speaker1, std_clean_speaker2, std_clean_speaker3, std_clean_speaker4];  % Error values
errors1 = [std(samp_mean1(:,1)), std(samp_mean1(:,2)), std(samp_mean1(:,3)), std(samp_mean1(:,4))];
errors2 = [std(samp_mean2(:,1)), std(samp_mean2(:,2)), std(samp_mean2(:,3)), std(samp_mean2(:,4))];
errors3 = [std(samp_mean3(:,1)), std(samp_mean3(:,2)), std(samp_mean3(:,3)), std(samp_mean3(:,4))];
errors4 = [std(samp_mean4(:,1)), std(samp_mean4(:,2)), std(samp_mean4(:,3)), std(samp_mean4(:,4))];

%Model Data
%clean_data = [];
% for m = 1:length(state.curfrgrid{end})
%     clean_data = [clean_data; state.curfrgrid{end}(1,:,m)];
% end


%errors2 = std(clean_data./max(mean(clean_data)));




% Create the plot with error bars
errorbar(x, y, errors, '-o','MarkerSize',8); hold on  % 'o' for markers at each data point
% plot(x,FRS1/max(FRS1));hold on
% plot(x,FRS2/max(FRS1));hold on
% plot(x,FRS3);hold on
% plot(x,FRS4);hold on



plot(x, y1, 'r-o','MarkerSize',8); hold on
%errorbar(x, y2, errors2, 'ro','MarkerSize',8); hold on
%errorbar(x, y3, errors3, 'go','MarkerSize',8); hold on
%errorbar(x, y4, errors4, 'ko','MarkerSize',8); hold on

%errorbar(x,mean(clean_data)./max(mean(clean_data)),errors2, 'o')
xlabel('X-axis label')  % Label for X-axis
ylabel('Y-axis label')  % Label for Y-axis
title('Plot with Error Bars')  % Title of the plot
grid on  % Add grid
xlim([0.5 4.5])

saveas(gcf,['FR' 'Masked' '.svg']);