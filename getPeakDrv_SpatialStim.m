% this script gets called by postProcessSims after spatial grid runs

% convert peakDrv to samples (10000 Hz)

close all;
load('peakDrv_target.mat');
figure(1);

% clearvars envelope
% fs = 44100; t = 1.5;
% envelope(:,1) = [ abs(sin(2*pi*AM_freqs(a)*(0:round(fs*1.5))'/fs)) ];
% envelope(:,2) = [ abs(sin(2*pi*AM_freqs(a)*(0:round(fs*1.5))'/fs)) ];
% stim_vec = (0:(size(envelope,1)-1))/fs;

t_bin = 20; %ms
t_vec = 0:t_bin:(3500-1);

pops = {'On','R1On','R2On'};
clearvars spks
for p = 1:length(pops)
    spks = sprintf('%s_V_spikes',pops{p});
    for a = 1:2
        figure(1);
        subplot(2,1,a);

        % convert peakDrv to ms and add 250ms (to account for zeropadding in Am
        % stimuli)
        peakDrv_ms = [round(peakDrv_target{a}(2:end),1)]+250;
        temp = [];
        % get spiketime indexes for each trial
        for t = (1:10)+(a-1)*10
            temp = [temp; find(snn_out(t).(spks))];
        end

        % convert to PSTH
        PSTH = histcounts(temp,0:t_bin*10:35000);
        plot(t_vec,PSTH,'linewidth',1); hold on;
        if p == length(pops)
            plot([peakDrv_ms;peakDrv_ms],[0 20]'.*ones(size(peakDrv_ms)),'--r','linewidth',1);
        end
        %     plot(stim_vec*1000,20*envelope(:,a),'b');
        ylabel(sprintf('Target %i Hz',a))

        % get activity within 150ms of each peakDrv event
        peakAct = [];
        for d = 1:length(peakDrv_ms)
            temp2 = temp/10 - peakDrv_ms(d);
            peak_vec = -150 : 6 : 150;
            peakAct(d,:) = histcounts(temp2,peak_vec);
        end
        figure(2);
        subplot(2,1,a);
        plot(peak_vec(1:end-1),mean(peakAct),'linewidth',1); hold on;
        if p ==  length(pops)
            plot([0 0],[0 10],'r');
        end
        ylabel(sprintf('Target %i Hz',a))
        %ylim([3 7])
    end
end

figure(1);
xlabel('Time (ms)');
legend(pops);
savefig(gcf,[simDataDir filesep 'AM response.fig']);

figure(2);
xlabel('Time from peakDrv (ms)');
legend(pops);
savefig(gcf,[simDataDir filesep 'peakDrv.fig']);

% %% Spike-triggered averages
% close all;
% 
% pops = {'R2On'};%{'On','R1On','R2On'};
% clearvars spks
% 
% fs_stim = 200000/5;
% y1_ds = downsample(y1,5)'; y2_ds = downsample(y2,5)';
% tlen_STA = .10; % [seconds]
% 
% figure;
% for p = 1:length(pops)
%     spks = sprintf('%s_V_spikes',pops{p});
%     for a = 1:2 % for each target
%         subplot(2,1,a);
%         temp = [];
%         % get spiketimes for each trial
%         for t = (1:10)+(a-1)*10
%             temp = [temp; find(snn_out(t).(spks))/10000];
%         end
% 
%         % for each spike time, get last 150 ms of target stimulus
% 
%         % convert to indexes from target
%         inds = round(temp * fs_stim); inds(inds < round(fs_stim*tlen_STA) | inds > length(y1_ds)) = [];
% 
%         t_STA = (-round(fs_stim*tlen_STA):0)/fs_stim;
%         inds = (-round(fs_stim*tlen_STA):0) + inds;
%         STA = [];
%         for i = 1:size(inds,1)
%             STA = cat(1,STA,eval(['y' num2str(a) '_ds(inds(i,:))']));
%         end
%         STA_mean = mean(STA);
% 
%         plot(t_STA,STA_mean,'linewidth',1); hold on;
%     end
% end
% sgtitle('Direct STA');
% 
% savefig(gcf,[simDataDir filesep 'direct STA.fig']);
% close;
% 
% %% STA from envelope
% 
% load('target_envelopes.mat','y_env');
% fs_stim = 200000/5;
% tlen_STA = .150; % [seconds]
% y1_env_ds = downsample(y_env{1}',5); y1_env_ds(y1_env_ds<0)=0;
% y2_env_ds = downsample(y_env{2}',5); y2_env_ds(y2_env_ds<0)=0;
% 
% y1_dB = 20*log10(y1_env_ds);
% y2_dB = 20*log10(y2_env_ds);
% 
% 
% pops = {'R2On'};%{'On','R1On','R2On'};
% clearvars spks
% 
% figure;
% for p = 1:length(pops)
%     spks = sprintf('%s_V_spikes',pops{p});
%     for a = 1:2 % for each target
%         subplot(2,1,a);
%         temp = [];
%         % get spiketimes for each trial
%         for t = (1:10)+(a-1)*10
%             temp = [temp; find(snn_out(t).(spks))/10000];
%         end
% 
%         % for each spike time, get last 150 ms of target stimulus
% 
%         % convert to indexes from target
%         inds = round(temp * fs_stim); inds(inds < round(fs_stim*tlen_STA) | inds > length(y1_env_ds)) = [];
% 
%         t_STA = (-round(fs_stim*tlen_STA):0)/fs_stim;
%         inds = (-round(fs_stim*tlen_STA):0) + inds;
%         STA = []; STA_dB=[];
%         for i = 1:size(inds,1)
%             STA = cat(1,STA,eval(['y' num2str(a) '_env_ds(inds(i,:))']));
%             STA_dB = cat(1,STA_dB,eval(['y' num2str(a) '_dB(inds(i,:))']));
%         end
%         STA_mean = mean(STA);
%         STA_dB(isinf(STA_dB)) = nan;
%         STA_dB_mean = mean(STA_dB,'omitnan');
% 
%         yyaxis left;
%         plot(t_STA,STA_mean,'linewidth',1);
%         yyaxis right;
%         plot(t_STA,STA_dB_mean,'linewidth',1);
%     end
% end
% sgtitle('Envelope STA');
% % legend('Linear','dB');
% % legend(pops)
% 
% savefig(gcf,[simDataDir filesep 'envelope STA.fig']);
