%Look accross all angles for the control.
%Isolate all ISIs
%Calculate max min.

%2. Isolate all ISIs in period of first onset
%Calculate max min


load('all_units_info_with_polished_criteria_modified_perf.mat','all_data');
load('sound_files.mat','sampleRate','target1','target2'); 

all_isi_whole_window = [];
all_isi_first_onset = [];

%4 angles
%220 cells
%10 trials
for k = 1:4
    for j = 1:220
        for m = 1:10
            %Do whole window
            if all_data(j).is_NS ~= 1
                spike_times_whole = all_data(j).ctrl_tar1_timestamps{m,k};
                mask = logical((spike_times_whole>0.1024) .* (spike_times_whole<0.2176));
                spike_times_window = spike_times_whole(mask);
                spike_times_whole_diff = diff(spike_times_whole);
                spike_times_window_diff = diff(spike_times_window);

                all_isi_whole_window = [all_isi_whole_window; spike_times_whole_diff];
                all_isi_first_onset = [all_isi_first_onset; spike_times_window_diff];

            end
        end
    end
end

%%

disp(['Minimum ISI for whole window across all (RS) cells for target 1 in clean condition: ' , num2str(min(all_isi_whole_window)*1000), ' ms'])
disp(['Maximum ISI for whole window across all (RS) cells for target 1 in clean condition: ' , num2str(max(all_isi_whole_window)*1000), ' ms'])
disp(['Minimum ISI for first onset across all (RS) cells for target 1 in clean condition: ' , num2str(min(all_isi_first_onset)*1000), ' ms'])
disp(['Maximum ISI for first onset across all (RS) cells for target 1 in clean condition: ' , num2str(max(all_isi_first_onset)*1000), ' ms'])

%%

figure(Position=[0,0,500,500]);
subplot(1,2,1)
histogram(all_isi_whole_window)
title('Whole window histogram')
ylabel('Frequency')
xlabel('ISI (s)')
xlim([0,1])
subplot(1,2,2)
histogram(all_isi_first_onset,'BinEdges',0:0.0001:0.1)
title('First onset window histogram')
ylabel('Frequency')
xlabel('ISI (s)')


%%

mean_whole = mean(all_isi_whole_window);
std_whole = std(all_isi_whole_window);
mean_windowed = mean(all_isi_first_onset);
std_windowed = std(all_isi_first_onset);
