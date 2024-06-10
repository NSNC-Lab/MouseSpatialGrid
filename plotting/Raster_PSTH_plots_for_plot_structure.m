close all
data_channel = 1;
fig_num = 6;

counter = 0;

cd Plot_Structure_Graphs\

titles = {'On','Off','ROn','ROff','SOnOff','TD','X','C'};


for j = 1:length(titles)
    plotting_data = [];
    for m = subz
        location = data(m).config;
        counter = counter + 1;
        for k = 1:10
            if k > 10
                d = 10;
            else
                d= 0;
            end

            %% Rasters
            %Grab the spikes that we are interested in
            on = transpose(data(m).spks.On(data_channel).channel1(k,:));
            off = transpose(data(m).spks.Off(data_channel).channel1(k,:));
            R1on = transpose(data(m).spks.ROn(data_channel).channel1(k,:));
            R2on = transpose(data(m).spks.ROff(data_channel).channel1(k,:));
            R1off = transpose(data(m).spks.SOnOff(data_channel).channel1(k,:));
            R2off = transpose(data(m).spks.TD(data_channel).channel1(k,:));
            S1OnOff = transpose(data(m).spks.X(data_channel).channel1(k,:));
            S2OnOff = transpose(data(m).spks.C(data_channel).channel1(k,:));
        
            reference_cell = {on,off,R1on,R2on,R1off,R2off,S1OnOff,S2OnOff};
            
            times = find(reference_cell{j} == 1);
            b = transpose(repmat(times,1,3));
            y_lines = nan(1,length(times));
            y_lines(1,:) = k-1-d;
            y_lines(2,:) = k-d;
            y_lines(3,:) = nan;
            figure(25+counter,'Visible', 'off')
            plotting_data = [plotting_data;times];
            plot(b,y_lines,'Color',[0 0 0 0.2],'LineWidth',0.3); hold on;

            

        end
        saveas(gcf, [titles{j} '_' location '_RASTER' '.fig']);

    

    
    raw_freq_data = histcounts(plotting_data, BinWidth=200);
    avg_data1 = movmean(raw_freq_data,3);

    %Conver to frequency from spike count
    % avg_data = (avg_data1*50/10)*10/150*75/40;
    
    figure(2000+counter,'Visible', 'off')
    plot(avg_data1,'r',LineWidth=0.5);

    saveas(gcf, ['titles{j}' '_' location '_PSTH' '.fig']); % Save as MATLAB figure file


    % if contains(titles{j}, 'S')
    %     plot(linspace(0,ending_sample-starting_sample,length(avg_data)),avg_data,'r',LineWidth=0.5);
    % else
    %     plot(linspace(0,ending_sample-starting_sample,length(avg_data)),avg_data,'k',LineWidth=0.5);
    % end
    % title(titles{j})
    % 
    % %Fig 4/5
    % xlim([3500 13500])
    % 
    % print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Desktop\Modeling Paper\Figures\Figure',num2str(fig_num),'\FpUpdate\',titles{j},'.svg']) % svg
    % 
    % 
    % %%Plot zoom raster & Voltage plot (Fig 5)
    % if fig_num == 5 || fig_num == 6
    % 
    %     indexs = find(snn_out(data_channel).R2On_V_spikes == 1);
    %     indexs2 = find(snn_out(data_channel).S2OnOff_V_spikes == 1);
    %     holder = snn_out(data_channel).R2On_V;
    %     holder(indexs) = 0;
    %     holder2 = snn_out(data_channel).S2OnOff_V;
    %     holder2(indexs2) = 0;
    % 
    % 
    %     if contains(titles{j}, 'S2OnOff')
    %         figure(90)
    %         plot(linspace(0,ending_sample-starting_sample,length(avg_data1)),avg_data1,'r',LineWidth=0.5); hold on
    %         %xlim([8500 13500])
    %         xlim([5500 10500])
    %         %print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Desktop\Modeling Paper\Figures\Figure',num2str(fig_num),'\FpUpdate\','ZoomedVersion','.svg']) % svg
    %         figure(91)
    %         %Range is changeable (Just looking for a relavent pattern here)
    %         plot(holder2(9000:15000),'r')
    %         ylim([-80 0])
    %         %print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Desktop\Modeling Paper\Figures\Figure',num2str(fig_num),'\FpUpdate\','PVspikes','.svg']) % svg
    % 
    %     elseif contains(titles{j}, 'R2on')
    %         figure(90)
    %         plot(linspace(0,ending_sample-starting_sample,length(avg_data1)),avg_data1,'k',LineWidth=0.5); hold on
    %         figure(92)
    %         plot(holder(9000:15000),'k')
    %         ylim([-80 0])
    %         %print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Desktop\Modeling Paper\Figures\Figure',num2str(fig_num),'\FpUpdate\','Espikes','.svg']) % svg
    %     end
    % 
    % 
    % end
    end

end

