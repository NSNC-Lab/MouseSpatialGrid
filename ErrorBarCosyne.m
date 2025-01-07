load('all_units_info_updated_manual.mat','all_data')

a = all_data(57).ctrl_tar1_masked_timestamps;
b = all_data(57).ctrl_tar1_timestamps;
nbins = 125;

%  this counts the number of spikes in 3 second, you need to devide it by 3
%  to get the FR
[count_ctrl_tar1_clean, ~, ~] = calc_num_spikes(b,nbins);
FR_ctrl_tar1_clean = count_ctrl_tar1_clean/3;

degree_list_x = [90 45 0 -90];
numBootstraps = 500;

% you need to modify the code "Calc_cen_mod_ERRF_bootstrap_raw"  by adding  the variable called "bootstrapMeans" as its output like this:
% function [errf_width_ctrl_all, mod_depth_ctrl_all, centroid_ctrl_all, bootstrapMeans]= calc_cen_mod_ERRF_bootstrap_raw(numBootstraps, counts,degree_list_x)  
[~, ~, ~, mean]= calc_cen_mod_ERRF_bootstrap_raw(numBootstraps, FR_ctrl_tar1_clean ,degree_list_x);  

% then you can then use "Calc_CI" to get the confidence interval at each speaker location.
CI_clean_speaker1 = calc_CI(mean(:,1), 0.95);
CI_clean_speaker2 = calc_CI(mean(:,2), 0.95);
CI_clean_speaker3 = calc_CI(mean(:,3), 0.95);
CI_clean_speaker4 = calc_CI(mean(:,4), 0.95);

% For masked condition, you need to know that "a = all_data(57).ctrl_tar1_masked_timestamps;"
% has a 10*4*4 structure, so the first 4 is target location and the second
% 4 is the masker location, they all follow the same order, 90 45 0 and -90

% so you could do this:
[count_ctrl_tar1_makser90, ~, ~] = calc_num_spikes(a(:,:,1),nbins);
FR_ctrl_tar1_masker90 = count_ctrl_tar1_makser90/3;