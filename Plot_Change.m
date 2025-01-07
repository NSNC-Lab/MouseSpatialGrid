

figure(54);

for var_change1 = [0:0.1:1]

    var_change2 = 0.2;
    var_change3 = 0.5;
    
    SpikingNetwork_paper;

    titles = {'on','off','R1on','R2on','R1Off','R2Off','S1OnOff','S2OnOff'};

    for j = 1:length(titles)
        for k = 1:10
            if k > 10
                d = 10;
            else
                d= 0;
            end
    
            %% Rasters
            %Grab the spikes that we are interested in
            on = transpose(data(subz).spks.On(data_channel).channel1(k,:));
            off = transpose(data(subz).spks.Off(data_channel).channel1(k,:));
            R1on = transpose(data(subz).spks.R1On(data_channel).channel1(k,:));
            R2on = transpose(data(subz).spks.R2On(data_channel).channel1(k,:));
            R1off = transpose(data(subz).spks.R1Off(data_channel).channel1(k,:));
            R2off = transpose(data(subz).spks.R2Off(data_channel).channel1(k,:));
            S1OnOff = transpose(data(subz).spks.S1OnOff(data_channel).channel1(k,:));
            S2OnOff = transpose(data(subz).spks.S2OnOff(data_channel).channel1(k,:));
    
            %reference_cell = {on,off,R1on,R2on,R1off,R2off,S1OnOff,S2OnOff};
            reference_cell = {on,off,R1on,R2on,R1off,R2off,S1OnOff,S2OnOff};
    
            
            
    
            times = find(reference_cell{j} == 1);
            b = transpose(repmat(times,1,3));
            y_lines = nan(1,length(times));
            y_lines(1,:) = k-1-d;
            y_lines(2,:) = k-d;
            y_lines(3,:) = nan;
            figure(25+j)
            plotting_data = [plotting_data;times];
            subplot(2,1,1);
            h2 = plot(b,y_lines,'Color',[0 0 0 0.3],'LineWidth',0.5); hold on
    
        end

    end

    
end