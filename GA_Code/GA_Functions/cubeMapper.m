% Example data
% p1 = linspace(0, 1, 100);
% p2 = linspace(0, 1, 100);
% p3 = linspace(0, 1, 100);
% 
% [Px, Py] = meshgrid(p1, p2);
% [Pz, Py2] = meshgrid(p3, p2);
% [Px2, Pz2] = meshgrid(p1, p3);
% 
% % Example heatmaps
% H1 = sin(Px) .* cos(Py); % Heatmap for (p1 & p2)
% H2 = sin(Px2) .* cos(Pz2); % Heatmap for (p1 & p3)
% H3 = sin(Pz) .* cos(Py2); % Heatmap for (p2 & p3)


%Set plot boundaries
y_bound = 0.3;
x_bound = 0.3;
Topo_plot_resolution = 100;
Num_Neighbors = 2;
lower_bound = 0;

var1 = 9;
var2 = 10;
var3 = 9;
var4 = 11;
var5 = 10;
var6 = 11;

%Track cumulative fitness at each area
total_fit12 = zeros(Topo_plot_resolution,Topo_plot_resolution);
total_fit13 = zeros(Topo_plot_resolution,Topo_plot_resolution);
total_fit23 = zeros(Topo_plot_resolution,Topo_plot_resolution);

% Animation loop
for gen = 1:length(state.curvars)

    %Clear the figures
    clf;
    hold on;
    axis([lower_bound x_bound lower_bound y_bound]);
    
    
    %General Plot
    %current_iteration =  reshape(state.curvars{gen},nVars,optionsGA.PopulationSize);
    current_iteration =  transpose([state.curvars{gen}(:,var1),state.curvars{gen}(:,var2)]);
    current_iteration2 =  transpose([state.curvars{gen}(:,var3),state.curvars{gen}(:,var4)]);
    current_iteration3 =  transpose([state.curvars{gen}(:,var5),state.curvars{gen}(:,var6)]);

    %Topo plot
    %1. Split the grid into a ton of little pieces.
    xEdges = linspace(lower_bound, x_bound, Topo_plot_resolution); 
    yEdges = linspace(lower_bound, y_bound, Topo_plot_resolution);

    all_fit = [];
    all_fit2 = [];
    all_fit3 = [];
    for i = 1:length(xEdges)
        columns_fit = [];
        columns_fit2 = [];
        columns_fit3 = [];
        for j = 1:length(yEdges)
            
            Distance_array = [];
            Distance_array2 = [];
            Distance_array3 = [];
            %Track the nearest nearest neightbors to each pixel on the
            %grid. Once nearest neighbors are found, take their avg.
            %fitness
            for k = 1:length(current_iteration)
                %L1
                %distFit = abs((current_iteration(1,k)-xEdges(i)))^2 + abs((current_iteration(2,k)-yEdges(j)))
                %L2
                distFit = sqrt(abs((current_iteration(1,k)-xEdges(i)))^2 + abs((current_iteration(2,k)-yEdges(j)))^2); 
                distFit2 = sqrt(abs((current_iteration2(1,k)-xEdges(i)))^2 + abs((current_iteration2(2,k)-yEdges(j)))^2); 
                distFit3 = sqrt(abs((current_iteration3(1,k)-xEdges(i)))^2 + abs((current_iteration3(2,k)-yEdges(j)))^2); 
                
                Distance_array = [Distance_array,distFit];
                Distance_array2 = [Distance_array2,distFit2];
                Distance_array3 = [Distance_array3,distFit3];
            end

            [~, sortedIndices] = sort(Distance_array);
            [~, sortedIndices2] = sort(Distance_array2);
            [~, sortedIndices3] = sort(Distance_array3);
            smallestIndices = sortedIndices(1:Num_Neighbors);
            smallestIndices2 = sortedIndices2(1:Num_Neighbors);
            smallestIndices3 = sortedIndices3(1:Num_Neighbors);
            
            accumulator = 0;
            accumulator2 = 0;
            accumulator3 = 0;
            for m = 1:Num_Neighbors
                accumulator = accumulator + state.curfitness{gen}(smallestIndices(m));
                accumulator2 = accumulator2 + state.curfitness{gen}(smallestIndices2(m));
                accumulator3 = accumulator3 + state.curfitness{gen}(smallestIndices3(m));
            end
            avg_fit = accumulator/Num_Neighbors;
            avg_fit2 = accumulator2/Num_Neighbors;
            avg_fit3 = accumulator3/Num_Neighbors;

            columns_fit = [columns_fit;avg_fit];
            columns_fit2 = [columns_fit2;avg_fit2];
            columns_fit3 = [columns_fit3;avg_fit3];

        end
        all_fit = [all_fit,columns_fit];
        all_fit2 = [all_fit2,columns_fit2];
        all_fit3 = [all_fit3,columns_fit3];

    end
    
    % Add a color bar to show the mapping
    colorbar;

    %pause(0.1);
    total_fit12 = total_fit12+all_fit;
    total_fit13 = total_fit13+all_fit2;
    total_fit23 = total_fit23+all_fit3;

end
hold off;
% Create figure
figure;


[Px, Py] = meshgrid(xEdges, yEdges);

% Plot (p1 & p2) on xy plane
surf(Px, Py, zeros(size(total_fit12)), total_fit12, 'EdgeColor', 'none');
hold on;

% Plot (p1 & p3) on xz plane
surf(Px, zeros(size(total_fit13)), Py, total_fit13, 'EdgeColor', 'none');

% Plot (p2 & p3) on yz plane
surf(zeros(size(total_fit23)), Px, Py, total_fit23, 'EdgeColor', 'none');

% Set view and labels
view(3);
xlabel('PV-sharpening');
ylabel('RC-boosting');
zlabel('SOM-Cross-channel');
colorbar;

% Set axis properties
axis tight;
grid on;
box on;
hold off;