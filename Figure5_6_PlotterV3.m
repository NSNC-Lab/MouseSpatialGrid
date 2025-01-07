close all

%5A : Bar chart showing the weights from On->pV and Off->PV %%%%%%%%%%%%%%%%%%%%%
figure(Position=[100,100,300,700]);

connection_type1 = 'On->S1OnOff';
connection_type2 = 'Off->S1OnOff';

con_idx1 = contains({varies.conxn},connection_type1) & contains({varies.param},'gSYN');
con_idx2 = contains({varies.conxn},connection_type2) & contains({varies.param},'gSYN');

on_PV_params = max(varies(con_idx1).range);  %Just using max here to compress this to 1 number, not actually taking a max here, although you could if you wanted
off_PV_params = max(varies(con_idx2).range);

values = [on_PV_params, off_PV_params];

barHandle = bar(values, 'FaceColor', 'none', 'EdgeColor', 'none');

% Hold the figure for overlaying lines
hold on;

x = 1:2;

% Manually draw lines for each bar
for i = x
    % X and Y coordinates for the bar edges
    x_left = x(i) - barHandle.BarWidth / 2; % Left edge
    x_right = x(i) + barHandle.BarWidth / 2; % Right edge
    y_top = values(i); % Top of the bar

    % Draw the left, right, and top edges
    plot([x_left, x_left], [0, y_top], 'k', 'LineWidth', 1.5); % Left edge
    plot([x_right, x_right], [0, y_top], 'k', 'LineWidth', 1.5); % Right edge
    plot([x_left, x_right], [y_top, y_top], 'k', 'LineWidth', 1.5); % Top edge
end


xticks(1:2);
xticklabels({'On->PV','Off->PV'});

set(gca, 'FontSize', 13, 'Box', 'off', 'LineWidth', 1);

xlabel('');
ylabel('g_{SYN} (\muS)', 'FontSize', 18);
xlim([0.3 2.7])
ylim([0 0.04]);
yticks(0:0.01:0.04);

print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(5),'\','5A','.svg'])


%5B : Raster and PSTH's for PV and E cells          %%%%%%%%%%%%%%%%%%


%Set how long you want the recording to be
starting_sample = 3500;
ending_sample = 10500;
%Fig5E_starting_sample = 3500;
%Fig5E_ending_sample = 13500;
num_iters = 30;

data_channel = 1;

titles = {'on','off','R1on','R2on','S1OnOff','S2OnOff'};


reg_grid = [1,2,3;4,5,6;7,8,9];
    
pos_grid_x = [0.1,0.1,0.1;0.4,0.4,0.4;0.7,0.7,0.7];
pos_grid_y = [0.7,0.45,0.2;0.7,0.45,0.2;0.7,0.45,0.2];

pos_grid_height = [0.2*ones(1,9)];
pos_grid_width = [0.27*ones(1,9)];

for j = 1:length(titles)
    

    if contains(titles{j},'on')
        col = 2;
    elseif contains(titles{j},'S')
        col = 1;
    else
        col = 3;
    end

    if contains(titles{j},'1')
        row = 2;
    elseif contains(titles{j},'2')
        row = 1;
    else 
        row = 3;
    end
   

    plotting_data = [];
    plotting_data2 = [];
    for k = 1:10
        if k > 10
            d = 10;
        else
            d= 0;
        end

        figure(10);
        
        al_idx = reg_grid(row,col);
        subplot('Position',[pos_grid_x(al_idx), pos_grid_y(al_idx), pos_grid_width(al_idx), pos_grid_height(al_idx)])

        %% Rasters
        %Grab the spikes that we are interested in
        on = transpose(data(subz).spks.On(data_channel).channel1(k,starting_sample:ending_sample));
        off = transpose(data(subz).spks.Off(data_channel).channel1(k,starting_sample:ending_sample));
        R1on = transpose(data(subz).spks.R1On(data_channel).channel1(k,starting_sample:ending_sample));
        R2on = transpose(data(subz).spks.R2On(data_channel).channel1(k,starting_sample:ending_sample));
        R1off = transpose(data(subz).spks.R1Off(data_channel).channel1(k,starting_sample:ending_sample));
        R2off = transpose(data(subz).spks.R2Off(data_channel).channel1(k,starting_sample:ending_sample));
        S1OnOff = transpose(data(subz).spks.S1OnOff(data_channel).channel1(k,starting_sample:ending_sample));
        S2OnOff = transpose(data(subz).spks.S2OnOff(data_channel).channel1(k,starting_sample:ending_sample));

        %reference_cell = {on,off,R1on,R2on,R1off,R2off,S1OnOff,S2OnOff};
        reference_cell = {on,off,R1on,R2on,S1OnOff,S2OnOff};

        times = find(reference_cell{j} == 1);
        b = transpose(repmat(times,1,3));
        y_lines = nan(1,length(times));
        y_lines(1,:) = k-1-d;
        y_lines(2,:) = k-d;
        y_lines(3,:) = nan;
        plotting_data = [plotting_data;times];
        h2 = plot(b,y_lines,'Color',[0 0 0 0.2],'LineWidth',0.5); hold on

        %TO DO Get AvgOf30 Activity
        %for data_channel = 1:30     
      
    end

    %TO DO plot the song underneath plot 8 or 9
    
    

    
    figure(10);

    set(gcf, 'Position', [100, 100, 800, 600]);

    subplot('Position',[pos_grid_x(al_idx), pos_grid_y(al_idx), pos_grid_width(al_idx), pos_grid_height(al_idx)])
    set(gca, 'Box', 'off'); 
    set(gca, 'XTick', []);  
    set(gca, 'YTick', []);

    if (contains(titles{j}, 'S2') || contains(titles{j}, 'R2'))
        raw_freq_data = histcounts(plotting_data, BinWidth=200,BinLimits=[0 ending_sample-starting_sample]);
    else
        raw_freq_data = histcounts(plotting_data, BinWidth=200,BinLimits=[0 ending_sample-starting_sample]);
    end

    
    avg_data1 = movmean(raw_freq_data,3);
    %avg_data = (avg_data1*50/10)*10/150;
    avg_data = (avg_data1*50/10)*10/150*75/40;
    if contains(titles{j}, 'S')
        %subplot(2,1,2);
        plot(linspace(0,ending_sample-starting_sample,length(avg_data)),avg_data,'r',LineWidth=1.5);
        %print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(4),'\',titles{j},'.svg'])
    else
        %subplot(2,1,2);
        plot(linspace(0,ending_sample-starting_sample,length(avg_data)),avg_data,'Color',[0.4660 0.6740 0.1880],LineWidth=1.5);
        %print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(4),'\',titles{j},'.svg'])
    end

    if contains(titles{j}, 'S2')
        PV_dat = avg_data;

    elseif contains(titles{j}, 'R2')
        E_dat = avg_data;
    end

    xlim([0 ending_sample-starting_sample])

subplot('Position',[0.7, 0.05, 0.27, 0.1])

[song1,fs] = audioread('200k_target1.wav');
plot(song1)
xlim([starting_sample*20 ending_sample*20])

set(gca, 'Box', 'off'); 
set(gca, 'XTick', []);  
set(gca, 'YTick', []);

subplot('Position',[0.4, 0.05, 0.27, 0.1])

%Song 1 operates 20x speed fo simulation 1/200k vs 1/10k for 0.1ms
%simulation step size
plot(song1)
xlim([starting_sample*20 ending_sample*20])
set(gca, 'Box', 'off'); 
set(gca, 'XTick', []);  
set(gca, 'YTick', []);


print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(5),'\','5B','.svg'])



end




%5C : Performance Bar Chart                                    %%%%%%%%%%%%%%%%%%%%%%%%
figure(Position=[300,100,500,700]);

values = [mean(perf.SPIKE), mean(perf.ISI), mean(perf.RISPIKE)];
errors = [std(perf.SPIKE), std(perf.ISI), std(perf.RISPIKE)];

barHandle = bar(values, 'FaceColor', 'none', 'EdgeColor', 'none');

% Hold the figure for overlaying lines
hold on;

x = 1:3;

% Manually draw lines for each bar
for i = x
    % X and Y coordinates for the bar edges
    x_left = x(i) - barHandle.BarWidth / 2; % Left edge
    x_right = x(i) + barHandle.BarWidth / 2; % Right edge
    y_top = values(i); % Top of the bar

    % Draw the left, right, and top edges
    plot([x_left, x_left], [0, y_top], 'k', 'LineWidth', 0.5); % Left edge
    plot([x_right, x_right], [0, y_top], 'k', 'LineWidth', 0.5); % Right edge
    plot([x_left, x_right], [y_top, y_top], 'k', 'LineWidth', 0.5); % Top edge
end

errorbar(x, values, errors, 'k.', 'LineWidth', 0.5, 'CapSize', 10);

xticks(1:3);
xticklabels({'SPIKE','ISI','RI-SPIKE'});

set(gca, 'FontSize', 13, 'Box', 'off', 'LineWidth', 1);

ax = gca; % Get current axes handle
ax.LineWidth = 0.5; % Set the axis line width
ax.TickLength = [0.02, 0.02]; % Adjust the tick size

xlabel('Spike distance measure', 'FontSize', 15);
ylabel('Performance', 'FontSize', 15);
xlim([0.3 3.7])
ylim([50 100]);
yticks(50:10:100);





print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(5),'\','5C','.svg'])


%5D : Voltage Plot            %%%%%%%%%%%%%%%%


%Need to bring back the voltage tracking in dynasim

indexs = find(snn_out(data_channel).R2On_V_spikes == 1);
indexs2 = find(snn_out(data_channel).S2OnOff_V_spikes == 1);
holder = snn_out(data_channel).R2On_V;
holder(indexs) = 0;
holder2 = snn_out(data_channel).S2OnOff_V;
holder2(indexs2) = 0;


figure;
subplot(2,1,1);
plot(holder2(starting_sample:ending_sample),'r', 'LineWidth', 1);
%xlim([0 2000])
subplot(2,1,2);
plot(holder(starting_sample:ending_sample),'k', 'LineWidth', 1);
%xlim([0 2000])

print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(5),'\','5D','.svg'])

%5E : PSTH Comparison Plot



%Representative
figure;
plot(linspace(0,ending_sample-starting_sample,length(avg_data)),PV_dat,'r',LineWidth=1); hold on
plot(linspace(0,ending_sample-starting_sample,length(avg_data)),E_dat,'k',LineWidth=1); hold on


print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(5),'\','5E','.svg'])
%Average
%figure;
%plot(linspace(0,ending_sample-starting_sample,length(avg_data)),PV_dat_avg,'r',LineWidth=0.5); hold on
%plot(linspace(0,ending_sample-starting_sample,length(avg_data)),E_dat_avg,'k',LineWidth=0.5); hold on



























%%%%ARCHIVE

        % %% Average of 30 rasters
        % on_avg = [];
        % off_avg = [];
        % R1on_avg = [];
        % R2on_avg = [];
        % R1off_avg = [];
        % R2off_avg = [];
        % S1OnOff_avg = [];
        % S2OnOff_avg = [];
        % 
        % for mk = 1:num_iters
        %     on_avg = [on_avg,transpose(data(subz).spks.On(mk).channel1(k,starting_sample:ending_sample))];
        %     off_avg = [off_avg,transpose(data(subz).spks.Off(mk).channel1(k,starting_sample:ending_sample))];
        %     R1on_avg = [R1on_avg,transpose(data(subz).spks.R1On(mk).channel1(k,starting_sample:ending_sample))];
        %     R2on_avg = [R2on_avg,transpose(data(subz).spks.R2On(mk).channel1(k,starting_sample:ending_sample))];
        %     R1off_avg = [R1off_avg,transpose(data(subz).spks.R1Off(mk).channel1(k,starting_sample:ending_sample))];
        %     R2off_avg = [R2off_avg,transpose(data(subz).spks.R2Off(mk).channel1(k,starting_sample:ending_sample))];
        %     S1OnOff_avg = [S1OnOff_avg,transpose(data(subz).spks.S1OnOff(mk).channel1(k,starting_sample:ending_sample))];
        %     S2OnOff_avg = [S2OnOff_avg,transpose(data(subz).spks.S2OnOff(mk).channel1(k,starting_sample:ending_sample))];
        % end
        % 
        % on_avg = mean(on_avg,2);
        % off_avg = mean(off_avg,2); 
        % R1on_avg = mean(R1on_avg,2); 
        % R2on_avg = mean(R2on_avg,2); 
        % R1off_avg = mean(R1off_avg,2); 
        % R2off_avg = mean(R2off_avg,2);
        % S1OnOff_avg = mean(S1OnOff_avg,2); 
        % S2OnOff_avg = mean(S2OnOff_avg,2);
        % 
        % on_holder = zeros(length(on_avg),1);
        % off_holder = zeros(length(off_avg),1);
        % R1on_holder = zeros(length(R1on_avg),1);
        % R2on_holder = zeros(length(R2on_avg),1);
        % R1off_holder = zeros(length(R1off_avg),1);
        % R2off_holder = zeros(length(R2off_avg),1);
        % S1OnOff_holder = zeros(length(R2off_avg),1);
        % S2OnOff_holder = zeros(length(S2OnOff_avg),1);
        % 
        % acc = 0;
        % q_start = 0;
        % for q = 1:length(on_avg)
        %     acc = acc + on_avg(q);
        %     if acc >= 1
        %         on_holder(round((q+q_start)/2)) = 1;
        %         q_start = q;
        %         acc = 0;
        %     end
        % end
        % 
        % 
        % acc = 0;
        % q_start = 0;
        % for q = 1:length(off_avg)
        %     acc = acc + off_avg(q);
        %     if acc >= 1
        %         off_holder(round((q+q_start)/2)) = 1;
        %         q_start = q;
        %         acc = 0;
        %     end
        % end
        % 
        % acc = 0;
        % q_start = 0;
        % for q = 1:length(R1on_avg)
        %     acc = acc + R1on_avg(q);
        %     if acc >= 1
        %         R1on_holder(round((q+q_start)/2)) = 1;
        %         q_start = q;
        %         acc = 0;
        %     end
        % end
        % 
        % acc = 0;
        % q_start = 0;
        % for q = 1:length(R2on_avg)
        %     acc = acc + R2on_avg(q);
        %     if acc >= 1
        %         R2on_holder(round((q+q_start)/2)) = 1;
        %         q_start = q;
        %         acc = 0;
        %     end
        % end
        % 
        % acc = 0;
        % q_start = 0;
        % for q = 1:length(R1off_avg)
        %     acc = acc + R1off_avg(q);
        %     if acc >= 1
        %         R1off_holder(round((q+q_start)/2)) = 1;
        %         q_start = q;
        %         acc = 0;
        %     end
        % end
        % 
        % acc = 0;
        % q_start = 0;
        % for q = 1:length(R2off_avg)
        %     acc = acc + R2off_avg(q);
        %     if acc >= 1
        %         R2off_holder(round((q+q_start)/2)) = 1;
        %         q_start = q;
        %         acc = 0;
        %     end
        % end
        % 
        % acc = 0;
        % q_start = 0;
        % for q = 1:length(S1OnOff_avg)
        %     acc = acc + S1OnOff_avg(q);
        %     if acc >= 1
        %         S1OnOff_holder(round((q+q_start)/2)) = 1;
        %         q_start = q;
        %         acc = 0;
        %     end
        % end
        % 
        % acc = 0;
        % q_start = 0;
        % for q = 1:length(S2OnOff_avg)
        %     acc = acc + S2OnOff_avg(q);
        %     if acc >= 1
        %         S2OnOff_holder(round((q+q_start)/2)) = 1;
        %         q_start = q;
        %         acc = 0;
        %     end
        % end
        % 
        % on_avg = on_holder;
        % off_avg = off_holder;
        % R1on_avg = R1on_holder;
        % R2on_avg = R2on_holder;
        % R1off_avg = R1off_holder;
        % R2off_avg = R2off_holder;
        % S1OnOff_avg = S1OnOff_holder;
        % S2OnOff_avg = S2OnOff_holder;
        % 
        % 
        % reference_cell = {on_avg,off_avg,R1on_avg,R2on_avg,R1off_avg,R2off_avg,S1OnOff_avg,S2OnOff_avg};
        % 
        % times = find(reference_cell{j} == 1);
        % b = transpose(repmat(times,1,3));
        % y_lines = nan(1,length(times));
        % y_lines(1,:) = k-1-d;
        % y_lines(2,:) = k-d;
        % y_lines(3,:) = nan;
        % plotting_data2 = [plotting_data2;times];
        % h2 = plot(b,y_lines,'Color',[0 0 0 0.3],'LineWidth',0.5); hold on


    %         figure(11);
    % subplot(3,3,reg_grid(row,col))
    % 
    % raw_freq_data = histcounts(plotting_data2, BinWidth=200,BinLimits=[0 ending_sample-starting_sample]);
    % avg_data1 = movmean(raw_freq_data,3);
    % %avg_data = (avg_data1*50/10)*10/150;
    % avg_data = (avg_data1*50/10)*10/150*75/40;
    % if contains(titles{j}, 'S')
    %     %subplot(2,1,2);
    %     plot(linspace(0,ending_sample-starting_sample,length(avg_data)),avg_data,'r',LineWidth=0.5);
    %     %print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(4),'\',titles{j},'.svg'])
    % else
    %     %subplot(2,1,2);
    %     plot(linspace(0,ending_sample-starting_sample,length(avg_data)),avg_data,'Color',[0.4660 0.6740 0.1880],LineWidth=0.5);
    %     %print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(4),'\',titles{j},'.svg'])
    % end
    % 
    % 
    % 
    % if contains(titles{j}, 'S2')
    %     PV_dat_avg = avg_data;
    % 
    % elseif contains(titles{j}, 'R2')
    %     E_dat_avg = avg_data;
    % end
    % 
    % title(titles{j})

