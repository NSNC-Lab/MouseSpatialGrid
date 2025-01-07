
figure;
%Set how long you want the recording to be
starting_sample = 1;
ending_sample = 35000;

data_channel = 1;
fig_num = 6;

%titles = {'on','off','R1on','R2on','R1Off','R2Off','S1OnOff','S2OnOff'};
titles = {'on'};
for j = 1:length(titles)
    plotting_data = [];
    for k = 1:120

        %% Rasters
        %Grab the spikes that we are interested in

        on = snn_out(k).On_V_spikes;
        off = snn_out(k).Off_V_spikes;
        R1on = snn_out(k).R1On_V_spikes;
        R2on = snn_out(k).R2On_V_spikes;
        R1off = snn_out(k).R1Off_V_spikes;
        R2off = snn_out(k).R2Off_V_spikes;
        S1OnOff = snn_out(k).S1OnOff_V_spikes;
        S2OnOff = snn_out(k).S2OnOff_V_spikes;

        %reference_cell = {on,off,R1on,R2on,R1off,R2off,S1OnOff,S2OnOff};
        %reference_cell = {on,off,R1on,R2on,R1off,R2off,S1OnOff,S2OnOff};
        reference_cell = {off};

        times = find(reference_cell{j} == 1);
        b = transpose(repmat(times,1,3));
        y_lines = nan(1,length(times));
        y_lines(1,:) = k-1;
        y_lines(2,:) = k;
        y_lines(3,:) = nan;
        figure(25+j)
        plotting_data = [plotting_data;times];
        subplot(4,1,1);
        h2 = plot(b,y_lines,'Color',[0 0 0 0.3],'LineWidth',0.5); hold on

    end


    xlim([1 length(off)])
    ylim([0 120])
    %ylim([0 20])
    
    raw_freq_data = histcounts(plotting_data, BinWidth=100,BinLimits=[0 ending_sample-starting_sample]);
    avg_data1 = movmean(raw_freq_data,3);
    %avg_data = (avg_data1*50/10)*10/150;
    avg_data = (avg_data1*50/10)*10/150*75/40;
    if contains(titles{j}, 'S')
        subplot(4,1,3);
        plot(linspace(0,ending_sample-starting_sample,length(avg_data)),avg_data,'r',LineWidth=0.5);
        %print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(4),'\',titles{j},'.svg'])
    else
        subplot(4,1,3);
        plot(linspace(0,ending_sample-starting_sample,length(avg_data)),avg_data,'Color',[0.4660 0.6740 0.1880],LineWidth=0.5);
        %print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(4),'\',titles{j},'.svg'])
    end
    title(titles{j})

    xlim([1 length(off)])

    subplot(4,1,2);

    [song1,~] = audioread('200k_target1.wav');

    %plot(song1)
    plot(full_stimuli)

    xlim([1 307200])
    %xlim([1 length(song1)])


    subplot(4,1,4)

    plot(fr_target_off{1,1})

    xlim([1 length(on)])

    
    %Fig 4/5
    %xlim([3500 13500])

    

    %print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Desktop\Modeling Paper\Figures\Figure',num2str(fig_num),'\FpUpdate\',titles{j},'.svg']) % svg


    

        %C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure4
    

end





% 
% 
% 
% 
%     end
% 
%     if k == 10
%         %ToDo insert the PSTH plot for each raster.
%     end
% 
% 
%     %a = transpose(data(subz).spks.R2Off(2).channel1(k,:));
%     %a2 = transpose(data(subz).spks.R1On(1).channel1(k,:));
%     %a = a(starting_sample:ending_sample);
%     %a2 = a2(starting_sample:ending_sample);
%     %times = find(a == 1);
%     %times2 = find(a2 == 1);
%     %b = transpose(repmat(times,1,3));
%     % b2 = transpose(repmat(times2,1,3));
%     % 
%     % if k > 10
%     %     d = 10;
%     % else
%     %     d= 0;
%     % end
% 
%     %plotting_data(k,:) = a;
%     %plotting_data2 = [plotting_data2;times];
%     %plotting_data3 = [plotting_data3;times2];
% 
%     %Start with prealocated size using nan
%     % y_lines = nan(1,length(times));
%     % y_lines(1,:) = k-1-d;
%     % y_lines(2,:) = k-d;
%     % y_lines(3,:) = nan;
%     % 
%     % y_lines2 = nan(1,length(times2));
%     % y_lines2(1,:) = k-1-d;
%     % y_lines2(2,:) = k-d;
%     % y_lines2(3,:) = nan;
%     % 
%     % 
%     % %Plot the rasters
%     % figure(25)
%     % h2 = plot(b,y_lines,'Color',[0 0 0 0.2],'LineWidth',0.3); hold on
%     % print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Desktop\Modeling Paper\Figures\Figure8\8HzRasterZoom','.svg']) % svg
%     % 
%     % xlim([1 35000])
%     %xlim([8500 20500])
%     %xlim([27000 30000])
%     %figure(13)
%     %h4 = plot(b2,y_lines2,'Color',[0 0 0 0.2],'LineWidth',0.3); hold on
%     %print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Desktop\Modeling Paper\Figures\Figure8\30HzRasterZoom','.svg']) % svg
%     %xlim([8500 20500])
%     %xlim([1 35000])
%     %xlim([27000 30000])
% 
% end
% 
% 
% %Adjust the size of the figure with figure handler (h)
% h.Position = [100 100 850 400];
% 
% %% PSTH
% %State if we are doing figure b or d
% b = true;
% %30 = 150hz
% 
% raw_freq_data = histcounts(plotting_data2, BinWidth=200,BinLimits=[0 ending_sample-starting_sample]);
% raw_freq_data2 = histcounts(plotting_data3, BinWidth=200,BinLimits=[0 ending_sample-starting_sample]);
% 
% avg_data = movmean(raw_freq_data,3);
% avg_data2 = movmean(raw_freq_data2,3);
% 
% figure(36);
% nbins = 175;
% histogram(plotting_data2,nbins,'FaceColor','k', 'FaceAlpha', 1);
% 
% figure(37);
% nbins = 175;
% histogram(plotting_data3,nbins,'FaceColor','k', 'FaceAlpha', 1);
% 
% % h = figure(25);
% % %hist_data = sum(plotting_data);
% % 
% % %histogram(plotting_data2, BinWidth=200)
% % if b == true
% %     %avg_data
% % 
% %     disp('h')
% %     %signal = movmean(sum(plotting_data),200);
% %     %sig_max = max(signal);
% %     %Scale
% %     avg_data = (avg_data*50/10)*10/150;
% %     %avg_data2 = (avg_data2*50/10)*10/150;
% % 
% % 
% % 
% % 
% %     %Create matching x-domain for PSTH
% %     xdom = linspace(0,ending_sample-starting_sample,length(avg_data));
% % 
% %     %plot psth
% %     %figure(99);
% %     plot(xdom,avg_data,'k',LineWidth=0.5); hold on
% %     %plot(xdom,avg_data2,'r',LineWidth=0.5);
% %     %xlim([8500 20500])
% %     xlim([3500 13500])
% %     %xlim([27000 30000])
% % else
% %     figure(2);
% %     plot(avg_data,'k',LineWidth=0.5)
% % end
