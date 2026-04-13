clear all
close all

%Load in the data and simulation
cd('C:\Users\ipboy\Desktop\SingleChannelFigures\JointMDS_Space')
load('50epoch_100ms_and_300ms_loss_All_cells.mat')
load('all_units_info_with_polished_criteria_modified_perf.mat','all_data');

%Extract out data and sim PSTHs
bin_width = 200;
time_len = 29801;
Rasters_sim = zeros(220,10,time_len);
PSTHs_sim = zeros([220,floor(time_len/bin_width)]);
[val,idx] = min(squeeze(losses(end,2,:,:))');

Rasters_data = zeros(220,10,time_len);
PSTHs_data = zeros([220,floor(time_len/bin_width)]);
bin_edges = 1:bin_width:time_len+1;
dummy_idxs = 1:time_len;

for k = 1:220
    Rasters_sim(k,:,:) = squeeze(output(k,idx(k),:,:,:));
    PSTHs_sim(k,:) = histcounts(squeeze(Rasters_sim(k,:,:)).*dummy_idxs,bin_edges);

    focus = 0;
    tuning = all_data(k).tuning_type;
    if strcmp(tuning, 'contra-tuned')
        focus = 1;
    elseif strcmp(tuning, '45°-tuned')
        focus = 2;
    elseif strcmp(tuning, 'center-tuned')
        focus = 3;
    elseif strcmp(tuning, 'ipsi-tuned')
        focus = 4;
    end

    SpikeTimes = all_data(k).ctrl_tar1_timestamps(:,focus);
    picture = zeros(10,29801);
    for m = 1:10
        stim_mask = logical((SpikeTimes{m} > 0) .* (SpikeTimes{m} < 2.9801));
        trial_indicies = round(SpikeTimes{m}(stim_mask)*10000);
        picture(m,trial_indicies+1) = 1;
    end
    Rasters_data(k,:,:) = picture;
    PSTHs_data(k,:) = histcounts(squeeze(Rasters_data(k,:,:)).*dummy_idxs,bin_edges);
end

%Joint MDS Analysis

%Distance Matrix calculation
comp_matrix = zeros([440,440]);
PSTHs_all = [PSTHs_data; PSTHs_sim];
for m = 1:440
    for k = 1:440
        comp_matrix(m,k) = sqrt(sum((PSTHs_all(m,:) - PSTHs_all(k,:)).^2));  %RMSE
    end
end

%Cells for labeling
cells = [1:220,1:220];

%MDS Plot (1)
[Y,stress] = mdscale(comp_matrix,2, 'criterion','metricsstress');
figure(Position=[0,0,1000,1000]);
scatter(Y(:,1),Y(:,2),36,'bo','filled'); hold on
labels = string(cells);
text(Y(:,1),Y(:,2),labels,'FontSize',10)
title('Full plot - Uncolored')
xlabel('MDS-Axis 1 (Firing Rate)')
ylabel('MDS-Axis 2 (Not Fully known)')

%MDS Plot (2) - Zoomed
figure(Position=[0,0,1000,1000]);
scatter(Y(:,1),Y(:,2),36,'bo','filled'); hold on
labels = string(cells);
text(Y(:,1),Y(:,2),labels,'FontSize',10)
title('Full plot - Uncolored')
xlabel('MDS-Axis 1 (Firing Rate)')
ylabel('MDS-Axis 2 (Not Fully known)')
ylim([-5 5])
xlim([-30 10])

%Color by self-similarity (normalized distance bewteen uper and lower
%matrix comparisons.
distances = [];
for k = 1:220
    dist_vals = 1./comp_matrix(220+k,1:220);
    max_dist = max(dist_vals);
    dist_vals_norm = dist_vals./max_dist;
    distances = [distances,dist_vals_norm(k)];
end
distances = [distances,distances];

%MDS Plot (3) - Colored
figure(Position=[0,0,1000,1000]);
scatter(Y(:,1),Y(:,2),36,distances,'filled'); hold on
labels = string(cells);
text(Y(:,1),Y(:,2),labels,'FontSize',10)
title('Full plot - Colored by distance Matrix')
xlabel('MDS-Axis 1 (Firing Rate)')
ylabel('MDS-Axis 2 (Not Fully known)')
colorbar;

%%

%Draw Lines between points
%MDS Plot (4) - Colored

f = figure(Position=[0,0,1000,1000]);
scatter(Y(:,1),Y(:,2),36,distances,'filled'); 
hold on

% Use whatever colormap you want
cmap = parula(256);

% If each line corresponds to one distance value, this should be length 220.
% Adjust this depending on how your distances are defined.
line_distances = distances(1:220);

% Set the color scaling
cmin = min(line_distances);
cmax = max(line_distances);
clim([cmin cmax])   % same color limits for scatter + lines

for k = 1:220
    % Normalize distance to [0,1]
    t = (line_distances(k) - cmin) / (cmax - cmin + eps);

    % Convert to colormap index
    idx = max(1, min(size(cmap,1), round(1 + t*(size(cmap,1)-1))));

    % RGB color for this line
    thisColor = cmap(idx,:);

    plot([Y(k,1),Y(k+220,1)], ...
         [Y(k,2),Y(k+220,2)], ...
         'Color', thisColor, 'LineWidth', 1.5);
end

labels = string(cells);
text(Y(:,1),Y(:,2),labels,'FontSize',10)

title('Full plot - Colored by distance Matrix')
xlabel('MDS-Axis 1 (Firing Rate)')
ylabel('MDS-Axis 2 (Not Fully known)')
colormap(cmap)
colorbar

exportgraphics(f, sprintf('Full_Plot_Colores.pdf'), 'ContentType','vector')

%%

%MDS Plot (5) - Zoomed

f = figure(Position=[0,0,1000,1000]);
scatter(Y(:,1),Y(:,2),36,distances,'filled'); 
hold on

cmap = parula(256);
line_distances = distances(1:220);

cmin = min(line_distances);
cmax = max(line_distances);
clim([cmin cmax])

for k = 1:220
    t = (line_distances(k) - cmin) / (cmax - cmin + eps);
    idx = max(1, min(size(cmap,1), round(1 + t*(size(cmap,1)-1))));
    thisColor = cmap(idx,:);

    plot([Y(k,1),Y(k+220,1)], ...
         [Y(k,2),Y(k+220,2)], ...
         'Color', thisColor, 'LineWidth', 1.5);
end

title('Full plot - Colored by distance Matrix')
xlabel('MDS-Axis 1 (Firing Rate)')
ylabel('MDS-Axis 2 (Not Fully known)')
colormap(cmap)
colorbar

% Set zoom limits BEFORE text
xl = [-35, -15];
yl = [-10, 10];
xlim(xl)
ylim(yl)

labels = string(cells);

% Only label points inside the visible window
mask = Y(:,1) >= xl(1) & Y(:,1) <= xl(2) & ...
       Y(:,2) >= yl(1) & Y(:,2) <= yl(2);

text(Y(mask,1), Y(mask,2), labels(mask), ...
    'FontSize', 10, 'Clipping', 'on')

exportgraphics(f, sprintf('Full_Plot_Colored_Zoomed.pdf'), ...
    'ContentType','vector')


%% Find verticality of the lines

verticality = [];

for k = 1:220
    verticality = [verticality,(180/pi)*atan(abs(Y(k,2)-Y(k+220,2))/abs(Y(k,1)-Y(k+220,1)))];
end

figure;
histogram(verticality)

disp('mean angle of lines: ' + string(mean(verticality)))


verticality = [];

boundary = 40;

for k = 1:220
    if sqrt(Y(k,1)^2 - Y(k,2)^2) > boundary
        if sqrt(Y(k+220,1)^2 - Y(k+220,2)^2) > boundary
            verticality = [verticality,(180/pi)*atan(abs(Y(k,2)-Y(k+220,2))/abs(Y(k,1)-Y(k+220,1)))];
        end
    end
end
figure;
histogram(verticality)

disp('mean angle of lines: ' + string(mean(verticality)))


%%

