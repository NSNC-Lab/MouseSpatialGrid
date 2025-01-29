%Select which plot
plot_num = 6;

%Plot the performance
figure(Position=[200,400,300,450]);
bar([mean(perf(1,:)),mean(perf(2,:)),mean(perf(3,:))],'FaceColor','none'); hold on
errorbar([1,2,3],[mean(perf(1,:)),mean(perf(2,:)),mean(perf(3,:))],[std(perf(1,:)),std(perf(1,:)),std(perf(1,:))],"LineStyle","none",'Color',[0,0,0]);
xlim([0.5,3.5])
ylim([50 100])
yticks([50:10:100])
xticklabels({'SPIKE','ISI','RI-SPIKE'})
ylabel('Performance')
xlabel('Spike distance measure')
box off;
print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(plot_num),'\Resubmission2025_2\','On_Conv_Performance_Timing_regime_High_PV','.svg']) % svg


%Plot the raster PSTH songs etc
starting_sample = 3500;
ending_sample = 13500;

figure(Position=[600,200,900,750]);
titles = {'on','off','S1OnOff','R1on'};
for j = 1:length(titles)
    plotting_data = [];

    subplot("Position",[0.15,j*0.17+0.11,0.35,0.15])

    %Song 1

    for k = 1:10
        if k > 10
            d = 10;
        else
            d= 0;
        end

        %% Rasters
        %Grab the spikes that we are interested in
        on = transpose(data(subz).spks.On(1).channel1(k,:));
        off = transpose(data(subz).spks.Off(1).channel1(k,:));
        S1OnOff = transpose(data(subz).spks.S1OnOff(1).channel1(k,:));
        R1on = transpose(data(subz).spks.R1On(1).channel1(k,:));
        
        reference_cell = {on,off,S1OnOff,R1on};

        times = find(reference_cell{j} == 1);
        b = transpose(repmat(times,1,3));
        y_lines = nan(1,length(times));
        y_lines(1,:) = k-1-d;
        y_lines(2,:) = k-d;
        y_lines(3,:) = nan;
        plotting_data = [plotting_data;times];
        h2 = plot(b,y_lines,'Color',[0 0 0 0.2],'LineWidth',0.3); hold on

    end
    
    raw_freq_data = histcounts(plotting_data, BinWidth=200,BinLimits=[starting_sample ending_sample]);
    avg_data1 = movmean(raw_freq_data,3);

    if contains(titles{j}, 'S')
        avg_data = avg_data1/10/0.02*(10/250); %Divide by 10 for # of trials, divide by 0.02 (bindWidth) to get spikes/s = Hz,next part is scale factor which will be reprented with scaling bar in paper
        plot(linspace(starting_sample,ending_sample,length(avg_data)),avg_data,'r',LineWidth=1.5); hold on
        yticklabels([0:125:250]);
        xticklabels('')
        
    else
        avg_data = avg_data1/10/0.02*(10/150); %Divide by 10 for # of trials, divide by 0.02 (bindWidth) to get spikes/s = Hz,next part is scale factor which will be reprented with scaling bar in paper
        plot(linspace(starting_sample,ending_sample,length(avg_data)),avg_data,'k',LineWidth=1.5); hold on
        yticklabels([0:75:150]);
        xticklabels('')
       
    end
    %title(titles{j})

    xlim([3500 13500])
    

    %Song 2

    plotting_data = [];

    subplot("Position",[0.55,j*0.17+0.11,0.35,0.15])

    for k = 11:20
        if k > 10
            d = 10;
        else
            d= 0;
        end

        %% Rasters
        %Grab the spikes that we are interested in
        on = transpose(data(subz).spks.On(1).channel1(k,:));
        off = transpose(data(subz).spks.Off(1).channel1(k,:));
        S1OnOff = transpose(data(subz).spks.S1OnOff(1).channel1(k,:));
        R1on = transpose(data(subz).spks.R1On(1).channel1(k,:));
        
        reference_cell = {on,off,S1OnOff,R1on};

        times = find(reference_cell{j} == 1);
        b = transpose(repmat(times,1,3));
        y_lines = nan(1,length(times));
        y_lines(1,:) = k-1-d;
        y_lines(2,:) = k-d;
        y_lines(3,:) = nan;
        plotting_data = [plotting_data;times];
        h2 = plot(b,y_lines,'Color',[0 0 0 0.2],'LineWidth',0.3); hold on

    end
    
    raw_freq_data = histcounts(plotting_data, BinWidth=200,BinLimits=[starting_sample ending_sample]);
    avg_data1 = movmean(raw_freq_data,3);

    if contains(titles{j}, 'S')
        avg_data = avg_data1/10/0.02*(10/250); %Divide by 10 for # of trials, divide by 0.02 (bindWidth) to get spikes/s = Hz,next part is scale factor which will be reprented with scaling bar in paper
        plot(linspace(starting_sample,ending_sample,length(avg_data)),avg_data,'r',LineWidth=1.5); hold on
        yticklabels([0:125:250]);
        xticklabels('')
        
        
    else
        avg_data = avg_data1/10/0.02*(10/150); %Divide by 10 for # of trials, divide by 0.02 (bindWidth) to get spikes/s = Hz,next part is scale factor which will be reprented with scaling bar in paper
        plot(linspace(starting_sample,ending_sample,length(avg_data)),avg_data,'k',LineWidth=1.5); hold on
        yticklabels([0:75:150]);
        xticklabels('')
        
        
    end
    %title(titles{j})

    xlim([3500 13500])
    box off;
    

end

subplot("Position",[0.15,0.15,0.35,0.1])
plot(song1(starting_sample*20:ending_sample*20),'k'); hold on
ylim([-1 1])
xlim([0 ending_sample*20-starting_sample*20])
yticklabels('')
xticks(0:(ending_sample*20-starting_sample*20)/5:ending_sample*20-starting_sample*20)
xticklabels(350:(1350-350)/5:1350)
xlabel('time (ms)')
box off;
subplot("Position",[0.55,0.15,0.35,0.1])
plot(song2(starting_sample*20:ending_sample*20),'k')
ylim([-1 1])
xlim([0 ending_sample*20-starting_sample*20])
yticklabels('')
xticks(0:(ending_sample*20-starting_sample*20)/5:ending_sample*20-starting_sample*20)
xticklabels(350:(1350-350)/5:1350)
xlabel('time (ms)')
box off;

text(-0.65, -0.7, 'Song 1', 'Units', 'normalized', ...
    'FontSize', 12, 'HorizontalAlignment', 'center');

text(0.5, -0.7, 'Song 2', 'Units', 'normalized', ...
    'FontSize', 12, 'HorizontalAlignment', 'center');

text(-1.4, 2, 'On', 'Units', 'normalized', ...
    'FontSize', 12, 'HorizontalAlignment', 'right');

text(-1.4, 3.75, 'Off', 'Units', 'normalized', ...
    'FontSize', 12, 'HorizontalAlignment', 'right');

text(-1.4, 5.5, 'PV1', 'Units', 'normalized', ...
    'FontSize', 12, 'HorizontalAlignment', 'right');

text(-1.4, 7.25, 'Output', 'Units', 'normalized', ...
    'FontSize', 12, 'HorizontalAlignment', 'right');




han = axes(gcf, 'visible', 'off'); % Create an invisible axes
han.YLabel.Visible = 'on';
ylabel(han, 'Spikes/s (Hz)');


print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(plot_num),'\Resubmission2025_2\','On_Conv_Rasters_Timing_regime_High_PV','.svg']) % svg


%Plot an average of 30 PSTHs between the two songs in order to show the
%discrimitability of the two songs.

titles = {'on','off','S1OnOff','R1on'};


song1_holder = [];
song2_holder = [];

E_ex = [];
I_ex = [];

starting_sample = 3500;
ending_sample = 8500;

for m = 1:length(FR)

    for j = 1:2
        plotting_data = [];
    
        %Song 1
    
        for k = 1:10
            if k > 10
                d = 10;
            else
                d= 0;
            end
    
            %% Rasters
            %Grab the spikes that we are interested in
            S1OnOff = transpose(data(subz).spks.S1OnOff(1).channel1(k,:));
            R1on = transpose(data(subz).spks.R1On(1).channel1(k,:));
            
            reference_cell = {S1OnOff,R1on};
    
            times = find(reference_cell{j} == 1);
            plotting_data = [plotting_data;times];
        end
        
        raw_freq_data = histcounts(plotting_data, BinWidth=200,BinLimits=[starting_sample ending_sample]);
        avg_data1 = movmean(raw_freq_data,3);
        avg_data = avg_data1/10/0.02*(10/250); %Divide by 10 for # of trials, divide by 0.02 (bindWidth) to get spikes/s = Hz,next part is scale factor which will be reprented with scaling bar in paper
        
        if j == 2
            song1_holder = [song1_holder;avg_data];
            E_ex = [E_ex;avg_data];
        else
            I_ex = [I_ex;avg_data];
        end
        
    
        %Song 2
    
        plotting_data = [];
    
        for k = 11:20
            if k > 10
                d = 10;
            else
                d= 0;
            end
    
            %% Rasters
            %Grab the spikes that we are interested in
            %Perhaps do the first layer. Some of the effects seem to be
            %deminished in the second set of PV cells. I think they are
            %necessary in order to match figure 4's depression rates.
            S1OnOff = transpose(data(subz).spks.S1OnOff(1).channel1(k,:));
            R1on = transpose(data(subz).spks.R1On(1).channel1(k,:));       
            reference_cell = {S1OnOff,R1on};
            times = find(reference_cell{j} == 1);
            plotting_data = [plotting_data;times];
    
        end
        
        raw_freq_data = histcounts(plotting_data, BinWidth=200,BinLimits=[starting_sample ending_sample]);
        avg_data1 = movmean(raw_freq_data,3);
        avg_data = avg_data1/10/0.02*(10/250); %Divide by 10 for # of trials, divide by 0.02 (bindWidth) to get spikes/s = Hz,next part is scale factor which will be reprented with scaling bar in paper
        
        if j == 2
            song2_holder = [song2_holder;avg_data];
        end
        
    
    end

end


figure;
subplot('Position',[0.1,0.5,0.8,0.4])
plot(mean(song1_holder,1),'b','LineWidth',1.5); hold on
plot(mean(song2_holder,1),'r','LineWidth',1.5); hold on
yticks(linspace(0,4,5))
yticklabels(linspace(0,100,5))
ylabel('Spikes/s (Hz)')
xticks('')

subplot('Position',[0.1,0.3,0.8,0.15])
plot(song1(starting_sample*20:ending_sample*20),'b'); hold on
xlim([0 ending_sample*20-starting_sample*20])
ylim([-1 1])
yticks('')
xticks('')
subplot('Position',[0.1,0.1,0.8,0.15])
plot(song2(starting_sample*20:ending_sample*20),'r'); hold on
xlim([0 ending_sample*20-starting_sample*20])
ylim([-1 1])
yticks('')
xticks(linspace(0,ending_sample*20-starting_sample*20,5))
xticklabels(linspace(350,850,5))
xlabel('time (ms)')

text(-0.02, 1.85, 'Song1', 'Units', 'normalized', ...
    'FontSize', 12, 'HorizontalAlignment', 'right');

text(-0.02, 0.5, 'Song2', 'Units', 'normalized', ...
    'FontSize', 12, 'HorizontalAlignment', 'right');

print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(plot_num),'\Resubmission2025_2\','On_Conv_Song1_Song2_Comp_Timing_regime_High_PV','.svg']) % svg


figure;
plot(mean(E_ex,1),'k','LineWidth',2); hold on
plot(mean(I_ex,1),'r','LineWidth',2); hold on
xticks(linspace(0,25,5))
xticklabels(linspace(350,850,5))
ylim([0 10])
yticks(linspace(0,10,3))
yticklabels(linspace(0,250,3))
xlabel('time (ms)')
ylabel('Spikes/s (Hz)')

print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(plot_num),'\Resubmission2025_2\','On_Conv_EI_Comparison_Timing_regime_High_PV','.svg']) % svg

