%For Average grid of a trial
% bestApproximateGrid2 = bestApproximateGrid;
%bestApproximateGrid = state.bestApproximateGridHistory{end};
% bestfrGrid2 = bestfrGrid;

bestfrGrid = mean(state.curfrgrid{2},3);
bestApproximateGrid = mean(state.curgrid{2},3);

% bestApproximateGrid = approximate_grid;
%bestfrGrid = state.frGridHistory{end};

% bestApproximateGrid = bestApproximateGrid2;
% 
% bestfrGrid = bestfrGrid2;


subz = 1:24;
locationLabels = {'90','45','0','-90'};

% setup populations
popNamesT = {'C'};
numPops = 1;        
popSizes = [1];
onlyC = true;

chanLabels = {'90','45','0','-90'};

% check if varied parameter is vector or matrix
numVars = 1;

for vv = 1:numVars
    %annotStr = annotTable{vv};

    % figure setup
    figwidth = max(length(popSizes) * 300, 900); % 300 * number of neurons?
    figheight = min(length(popSizes) * 300, 800);
    h = figure('position', [200 50 figwidth*1.5 figheight*2]);

    plotheight = 1 / numPops * 0.55; % percentage
    xstart = 0.1;
    ystart = 0.1;


    for pop = 1:numPops
        % vary grid width depending on the number of neurons in subpopulation
        if popSizes(pop) == 1
            plotwidth = 0.233; % percentage, 0.7*num neurons
        else
            plotwidth = 0.14;
        end

        %for chan = 1:popSizes(pop)
        for chan = 1:1 %Changed 7/19 to just look for upper left corner IB
            
            % subfigure positions
            col = chan - 1;
            row = pop - 1;
            xoffset = xstart + plotwidth * (col) * 1.2;
            yoffset = ystart + plotheight * 1.6 * row;

            % mixed cases

                for i = 1:1
                    
                    perf.C = bestApproximateGrid(2:5,:);
                    fr.C= bestfrGrid(2:5,:);
                end
                ax = subplot('Position', [xoffset yoffset plotwidth plotheight]);
                plotPerfGrid(flipud(perf.C), flipud(fr.C), []);
                

                % add axes
                if onlyC || (chan == 1 && pop == 1)
                    xticks(1:4);
                    xticklabels(locationLabels);
                    yticks(1:4);
                    yticklabels(fliplr(locationLabels));
                    xlabel('target location');
                    ylabel('masker location');
                end


            % target-only configs
           

    
    
                perf.C2 = bestApproximateGrid(1,:);
                fr.C2 = bestfrGrid(1,:);

                ax = subplot('Position', [xoffset yoffset + plotheight + 0.01 plotwidth plotheight * 0.2]);



                plotPerfGrid(perf.C2, fr.C2, '');
                

              
        end
    end

end

figure(1000)
plot(bestfrGrid(1,:)); hold on; %Clean\
plot(bestfrGrid(2,:)); hold on; %Mixed Masker at -90


%Target
plot(Target_fr_grid(1,:),'--'); hold on
plot(Target_fr_grid(2,:),'--')


legend({'Clean Approx','Mixed Approx','Clean Target','Mixed Target'})