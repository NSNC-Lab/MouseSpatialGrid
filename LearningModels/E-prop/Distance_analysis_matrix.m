%% 1) For all of the outputs extract a PSTH
cd('C:\Users\ipboy\Documents\GitHub\ModelingEffort\Multi-Channel\Plotting\OliverDataPlotting')
load('output_compressed_Eprop_All_cells_20260328_203631')
load('all_units_info_with_polished_criteria_modified_perf.mat','all_data');
[target1, fs] = audioread('C:\Users\ipboy\Documents\GitHub\MouseSpatialGrid\LearningModels\resampled-stimuli\target\200k_target1.wav');

%a Extract the best fit from all of the cells

[a,idx] = min(squeeze(losses(100,2,:,:))'); %Idx is the index of the best amoung the batch

%Look at an example
figure(position=[0,0,1000,1000]);%Data plot
spy(squeeze(output(105,idx(105),:,:,:)))

rasters = zeros([220,10,29801]);
for k = 1:220
    rasters(k,:,:,:) = squeeze(output(k,idx(k),:,:,:));
end

bin_width = 200; %In 0.1 ms
bin_edges = 1:bin_width:size(rasters,3)+1;

PSTHs = zeros([220,149]);
dummy_idxs = 1:size(rasters,3);
for k = 1:220
    PSTHs(k,:) = histcounts(squeeze(rasters(k,:,:)).*dummy_idxs,bin_edges);
end

matrix = zeros([220,220]);

for m = 1:220
    for k = 1:220
        [r,lags] = xcorr(PSTHs(m,:), PSTHs(k,:), 'coeff');
        matrix(m,k) = max(r);
    end
end

figure;
h = pcolor(matrix);  %Looks like pcolor is flipping things to make it like a spectorgram (Just an fyi)
set(h, 'edgecolor', 'none')

%% 2) For all of the data extract a PSTH

data_rasters = zeros([220,10,29801]);
for n = 1:220
    
    SpikeTimes = all_data(n).ctrl_tar1_timestamps(:,1);
    picture = zeros(10,29801);
    
    for m = 1:10
        stim_mask = logical((SpikeTimes{m} > 0) .* (SpikeTimes{m} < 2.9801));
        trial_indicies = round(SpikeTimes{m}(stim_mask)*10000);
        picture(m,trial_indicies+1) = 1;
    end

    data_rasters(n,:,:) = picture;
end

%Look at an example
figure(position=[0,0,1000,1000]);%Data plot
spy(squeeze(data_rasters(105,:,:)))

data_PSTHs = zeros([220,149]);
for k = 1:220
    data_PSTHs(k,:) = histcounts(squeeze(data_rasters(k,:,:)).*dummy_idxs,bin_edges);
end

data_matrix = zeros([220,220]);

for m = 1:220
    for k = 1:220
        [r,lags] = xcorr(data_PSTHs(m,:), data_PSTHs(k,:), 'coeff');
        data_matrix(m,k) = max(r);
    end
end

figure;
h = pcolor(data_matrix);  %Looks like pcolor is flipping things to make it like a spectorgram (Just an fyi)
set(h, 'edgecolor', 'none')



%% 3) Combine the two PSTH holders

PSTH_all = [PSTHs;data_PSTHs];
all_matrix = zeros([220,220]);

for m = 1:440
    for k = 1:440
        [r,lags] = xcorr(PSTH_all(m,:), PSTH_all(k,:), 'coeff');
        all_matrix(m,k) = max(r);
    end
end

%% 4) Color code this matrix

figure;
h = pcolor(all_matrix);  %Looks like pcolor is flipping things to make it like a spectorgram (Just an fyi)
set(h, 'edgecolor', 'none')


%% 5) Look at data vs sim and normalize by row max.
comp_matrix = zeros([220,220]);
for m = 1:220
    for k = 1:220
        %[r,lags] = xcorr(data_PSTHs(m,:), PSTHs(k,:), 'coeff');
        %r = fitlm(data_PSTHs(m,:), PSTHs(k,:)).Rsquared.Ordinary;
        r = sum((data_PSTHs(m,:) - PSTHs(k,:)).^2);
        %comp_matrix(m,k) = max(r);
        comp_matrix(m,k) = r;
    end
    max_val = max(comp_matrix(m,:));
    comp_matrix(m,:) = 1./(comp_matrix(m,:)/max_val);
    max_val = max(comp_matrix(m,:));
    comp_matrix(m,:) = (comp_matrix(m,:)/max_val);
end


figure;
h = pcolor(comp_matrix);  %Looks like pcolor is flipping things to make it like a spectorgram (Just an fyi)
set(h, 'edgecolor', 'none')


%% 6) MSE shows up in sim vs data. Plot whole matrix 

%It doesn't show up well in the full matrix. Gets drowned out by the center
%line.

%Try without the centerline
%Line pops without cernterline

comp_matrix2 = zeros([440,440]);

new_mat = zeros([440,439]);

for m = 1:440
    for k = 1:440
        %if m ~= k
        %[r,lags] = xcorr(data_PSTHs(m,:), PSTHs(k,:), 'coeff');
        %r = fitlm(data_PSTHs(m,:), PSTHs(k,:)).Rsquared.Ordinary;
        r = sum((PSTH_all(m,:) - PSTH_all(k,:)).^2);
        %comp_matrix(m,k) = max(r);
        comp_matrix2(m,k) = r; %NOTE! Adding an offset just for visualization purposes only... Other wise we divide by 0
        %end
    end
    
    %if m == 118
    %    disp('')
    %end
    new_mat(m,:) = comp_matrix2(m,m~=(1:440));

    max_val = max(new_mat(m,:));
    new_mat(m,:) = 1./(new_mat(m,:)/max_val);
    max_val = max(new_mat(m,:));
    new_mat(m,:) = (new_mat(m,:)/max_val);
end




figure;
h = pcolor(new_mat);  %Looks like pcolor is flipping things to make it like a spectorgram (Just an fyi)
set(h, 'edgecolor', 'none')

%% Giant Spy plot?
figure;
subplot(1,2,1);
spy(reshape(Rasters_100_300ms,220*10,29801))
axis normal
title('sim')
subplot(1,2,2);
spy(reshape(Rasters_data,220*10,29801))
axis normal
title('data')


%% Lets looks at MDS plot

%Best Criterions
%1. metricsstress
%2. Sammon
%3. metricstresss


comp_matrix = zeros([220,220]);
for m = 1:220
    for k = 1:220
        %[r,lags] = xcorr(data_PSTHs(m,:), PSTHs(k,:), 'coeff');
        %r = fitlm(data_PSTHs(m,:), PSTHs(k,:)).Rsquared.Ordinary;
        r = sum((data_PSTHs(m,:) - data_PSTHs(k,:)).^2);
        %comp_matrix(m,k) = max(r);
        comp_matrix(m,k) = r;
    end
    %max_val = max(comp_matrix(m,:));
    %comp_matrix(m,:) = 1./(comp_matrix(m,:)/max_val);
    %max_val = max(comp_matrix(m,:));
    %comp_matrix(m,:) = (comp_matrix(m,:)/max_val);
end

%Y1 = mdscale(comp_matrix,2, 'criterion','sammon');
%Y2 = mdscale(comp_matrix,2, 'criterion','strain');
[Y3,stress] = mdscale(comp_matrix,2, 'criterion','metricsstress');

%Use shepards plot to assess goodnes of fit? According to matlab mdscale
%documentation
figure;
%distances1 = pdist(Y1);
%distances2 = pdist(Y2);
distances3 = pdist(Y3);
upper_no_diag = comp_matrix(tril(true(size(comp_matrix)),-1));
%plot(upper_no_diag,distances1,'bo', [0 200],[0 200],'k--'); hold on
%plot(upper_no_diag,distances2,'rx', [0 200],[0 200],'k--'); hold on
plot(upper_no_diag,distances3,'bo', [0 200],[0 200],'k--'); hold on
plot(1:25000,1:25000,'k--','LineWidth',2)
xlabel('Dissimilarities')
ylabel('Distances')


figure;
plot(Y3(:,1),Y3(:,2),'ko');

%% Quantify the onsetness of the data

target_clipped = target1(50000:end);

[yupper,ylower] = envelope(target_clipped,5000,'peak');
figure;
plot(target_clipped); hold on
plot(yupper)
xlim([18500 582000])

figure;
subplot(4,1,1)
plot(yupper>0.1)
ylim([0,2])
subplot(4,1,2)
plot(yupper>0.13)
ylim([0,2])
subplot(4,1,3)
% 0.16 looks like it captures the main onsets pretty well
plot(yupper>0.16)
ylim([0,2])
subplot(4,1,4)
plot(yupper>0.2)
ylim([0,2])


%Lets use 0.16 and quantify onsetness as FR_during_onsets/FR_during_offsets
mask = (yupper>0.16) .* (1:(646011-50000+1) > 18000)'; %Remove the artifact at the begining

%figure;
%plot(mask)

%Donwsample mask to the PSTH timescale (149 bins)
mask_PSTH = downsample(mask,ceil((646011-50000+1)/149));
mask_PSTH_numbered = (1:149).*mask_PSTH';

mask_PSTH_ordinal = mask_PSTH_numbered(((1:149).*mask_PSTH_numbered)>0);

%Offset
mask_PSTH_offset = 1-mask_PSTH;
mask_PSTH_numbered_offset = (1:149).*mask_PSTH_offset';
mask_PSTH_ordinal_offset = mask_PSTH_numbered_offset(((1:149).*mask_PSTH_numbered_offset)>0);

%figure;
%plot(mask_PSTH)

onsetness = sum(data_PSTHs(:,mask_PSTH_ordinal),2)./sum(data_PSTHs(:,mask_PSTH_ordinal_offset),2);

%We'll do non-strict asignments per PSTH bin


figure;
spy(squeeze(data_rasters(203,:,:)))


%color by onsetnes
figure;
scatter(Y3(:,1),Y3(:,2),10,layer,'filled');
colorbar

%% look at the same analysis for the fits

comp_matrix = zeros([220,220]);
for m = 1:220
    for k = 1:220
        %[r,lags] = xcorr(data_PSTHs(m,:), PSTHs(k,:), 'coeff');
        %r = fitlm(data_PSTHs(m,:), PSTHs(k,:)).Rsquared.Ordinary;
        r = sum((PSTHs(m,:) - PSTHs(k,:)).^2);
        %comp_matrix(m,k) = max(r);
        comp_matrix(m,k) = r;
    end
    %max_val = max(comp_matrix(m,:));
    %comp_matrix(m,:) = 1./(comp_matrix(m,:)/max_val);
    %max_val = max(comp_matrix(m,:));
    %comp_matrix(m,:) = (comp_matrix(m,:)/max_val);
end

%Y1 = mdscale(comp_matrix,2, 'criterion','sammon');
%Y2 = mdscale(comp_matrix,2, 'criterion','strain');
[Y3,stress] = mdscale(comp_matrix,2, 'criterion','metricsstress');
Y4 = Y3;

%Use shepards plot to assess goodnes of fit? According to matlab mdscale
%documentation
figure;
%distances1 = pdist(Y1);
%distances2 = pdist(Y2);
distances3 = pdist(Y3);
upper_no_diag = comp_matrix(tril(true(size(comp_matrix)),-1));
%plot(upper_no_diag,distances1,'bo', [0 200],[0 200],'k--'); hold on
%plot(upper_no_diag,distances2,'rx', [0 200],[0 200],'k--'); hold on
plot(upper_no_diag,distances3,'bo', [0 200],[0 200],'k--'); hold on
plot(1:25000,1:25000,'k--','LineWidth',2)
xlabel('Dissimilarities')
ylabel('Distances')


figure;
plot(Y3(:,1),Y3(:,2),'ko');

%figure;
%plot(mask)


%figure;
%plot(mask_PSTH)

onsetness = sum(data_PSTHs(:,mask_PSTH_ordinal),2)./sum(data_PSTHs(:,mask_PSTH_ordinal_offset),2);

layer = []

for k = 1:220
    if  strcmp(all_data(k).layer,'L2/3')
        layer = [layer,1];
    elseif strcmp(all_data(k).layer, 'L5/6')
        layer = [layer,2];
    elseif strcmp(all_data(k).layer,'L4')
        layer = [layer,3];
    elseif strcmp(all_data(k).layer, 'NaN')
        layer = [layer,4];
    end
end

%We'll do non-strict asignments per PSTH bin


%figure;
%spy(squeeze(data_rasters(203,:,:)))


%color by onsetnes
figure;
scatter(Y3(:,1),Y3(:,2),10,layer,'filled');
colorbar


%% Procruste
[d,Z,transform] = procrustes(Y3,Y4);


%% Save out for futher anlaysis

save("analysis_data.mat","PSTHs","data_PSTHs","mask_PSTH_ordinal")
save("analysis_data2.mat","rasters","data_rasters")


%% Kamal Suggestion (latest)
% Are you currently using the MSE as a measure of distance? It might be better to use the square root of MSE for metric MDS: 
% this is more akin to Euclidean distance.
%Let's try generating the MDS plot for the data and the model in separate figures (as you did) using two different coloring schemes.
%a) every cell is a different color
%            b) color each cell by layer (2/3, 4 and 5/6).
%  3. Try plotting the figures on a log log axis (I think it's easy to change this in MATLAB). 
% This may make it easier to visualize clusters in dense regions of the plot and still allow us to see the full range of the dat


%Data MDS plot

comp_matrix_data = zeros([220,220]);
comp_matrix_sim = zeros([220,220]);
for m = 1:220
    for k = 1:220
        comp_matrix_data(m,k) = sqrt(sum((data_PSTHs(m,:) - data_PSTHs(k,:)).^2));  %RMSE
        comp_matrix_sim(m,k) = sqrt(sum((PSTHs(m,:) - PSTHs(k,:)).^2));  
    end
end

[Y_data,stress_data] = mdscale(comp_matrix_data,2, 'criterion','metricsstress');
[Y_sim,stress_sim] = mdscale(comp_matrix_sim,2, 'criterion','metricsstress');

layer = [];

for k = 1:220
    if  strcmp(all_data(k).layer,'L2/3')
        layer = [layer,1];
    elseif strcmp(all_data(k).layer, 'L5/6')
        layer = [layer,2];
    elseif strcmp(all_data(k).layer,'L4')
        layer = [layer,3];
    elseif strcmp(all_data(k).layer, 'NaN')
        layer = [layer,4];
    end
end

figure(position=[0,0,1000,1000]);%Data plot
xplot = sign(Y_data(:,1)).*log10(1+abs(Y_data(:,1)));
yplot = sign(Y_data(:,2)).*log10(1+abs(Y_data(:,2)));
scatter(xplot,yplot,10,layer); hold on
labels = string(1:220);
text(xplot,yplot,labels,'FontSize',10)
%set(gca,'XScale','log','YScale','log') %Negetive values. Can still do log
%compression but will need to finagle it.
xlabel('mds axis 1')
ylabel('mds axis 2')

figure(position=[0,0,1000,1000]); %Sim plot
xplot = sign(Y_sim(:,1)).*log10(1+abs(Y_sim(:,1)));
yplot = sign(Y_sim(:,2)).*log10(1+abs(Y_sim(:,2)));
scatter(xplot,yplot,10,layer); hold on
labels = string(1:220);
text(xplot,yplot,labels,'FontSize',10)
xlabel('mds axis 1')
ylabel('mds axis 2')

figure(position=[0,0,1000,1000]); %Sim plot (y-axis flipped)
xplot = sign(Y_sim(:,1)).*log10(1+abs(Y_sim(:,1)));
yplot = sign(Y_sim(:,2)).*log10(1+abs(Y_sim(:,2)));
scatter(xplot,-yplot,10,layer); hold on
labels = string(1:220);
text(xplot,-yplot,labels,'FontSize',10)
xlabel('mds axis 1')
ylabel('mds axis 2')

%Pull from each "cluster"
%Upper right cluster
figure;
subplot(4,1,1)
spy(squeeze(data_rasters(136,:,:)))
title('data - cell 136')
subplot(4,1,2)
spy(squeeze(rasters(136,:,:)))
title('sim - cell 136')
subplot(4,1,3)
plot(data_PSTHs(136,:)); hold on
plot(PSTHs(136,:)); 
legend({'data','sim'})
title('PSTH overlay')
subplot(4,1,4)
plot(target1)
title('stimulus')

%Dense center cluster
cell = 23;
figure;
subplot(4,1,1)
spy(squeeze(data_rasters(cell,:,:)))
title('data - cell ' + string(cell))
subplot(4,1,2)
spy(squeeze(rasters(cell,:,:)))
title('sim - cell ' + string(cell))
subplot(4,1,3)
plot(data_PSTHs(cell,:)); hold on
plot(PSTHs(cell,:)); 
legend({'data','sim'})
title('PSTH overlay')
subplot(4,1,4)
plot(target1)
title('stimulus')

%Lower Right Cluster
cell = 202;
figure(position=[0,0,1000,1000]);
subplot(4,1,1)
spy(squeeze(data_rasters(cell,:,:)))
title('data - cell ' + string(cell))
subplot(4,1,2)
spy(squeeze(rasters(cell,:,:)))
title('sim - cell ' + string(cell))
subplot(4,1,3)
plot(data_PSTHs(cell,:)); hold on
plot(PSTHs(cell,:)); 
legend({'data','sim'})
title('PSTH overlay')
subplot(4,1,4)
plot(target1(50000:end))
title('stimulus')

%% Color the MDS plots by onsetness

%look at skewed onsetness to see extremes better (and unity)

skew_onsetness = onsetness.^(1/2);

figure(position=[0,0,1000,1000]); %Data plot
xplot = sign(Y_data(:,1)).*log10(1+abs(Y_data(:,1)));
yplot = sign(Y_data(:,2)).*log10(1+abs(Y_data(:,2)));
scatter(xplot,yplot,72,layer,'filled'); hold on
labels = string(1:220);
%text(xplot,yplot,labels,'FontSize',10)
%set(gca,'XScale','log','YScale','log') %Negetive values. Can still do log
%compression but will need to finagle it.
xlabel('mds axis 1')
ylabel('mds axis 2')
colorbar

figure(position=[0,0,1000,1000]); %Sim plot (y-axis flipped)
xplot = sign(Y_sim(:,1)).*log10(1+abs(Y_sim(:,1)));
yplot = sign(Y_sim(:,2)).*log10(1+abs(Y_sim(:,2)));
scatter(xplot,-yplot,72,layer,'filled'); hold on
labels = string(1:220);
%text(xplot,-yplot,labels,'FontSize',10)
xlabel('mds axis 1')
ylabel('mds axis 2')
colorbar

%% Comapre 200ms with 100ms  3.30.2026
clear all
close all
load('all_units_info_with_polished_criteria_modified_perf.mat','all_data');
load('output_compressed_Eprop_All_cells_20260328_203631')

losses_100ms = losses;
output_100ms = output;
params_100ms = params;

load('output_compressed_Eprop_All_cells_20260324_042421')

losses_200ms = losses;
output_200ms = output;
params_200ms = params;

load('output_compressed_Eprop_All_cells_20260331_060643')

losses_100_300ms = losses;
output_100_300ms = output;
params_100_300ms = params;

clear losses output params 

% Compare the losses
figure;
subplot(9,1,1)
plot(squeeze(mean(losses_100ms(:,2,:,:),4)))
title('100ms loss full view')
subplot(9,1,2)
plot(squeeze(mean(losses_100ms(:,2,:,:),4)))
ylim([0, 40000])
title('100ms zoomed')
subplot(9,1,3)
plot(squeeze(mean(mean(losses_100ms(:,2,:,:),4),3)))
ylim([0, 40000])
title('100ms mean_all')
subplot(9,1,4)
plot(squeeze(mean(losses_200ms(:,2,:,:),4)))
xlim([0, 100])
title('200ms loss full view')
subplot(9,1,5)
plot(squeeze(mean(losses_200ms(:,2,:,:),4)))
ylim([0, 40000])
xlim([0, 100])
title('200ms zoomed')
subplot(9,1,6)
plot(squeeze(mean(mean(losses_200ms(:,2,:,:),4),3)))
ylim([0, 40000])
xlim([0, 100])
title('200ms mean_all')
subplot(9,1,7)
plot(squeeze(mean(losses_100_300ms(:,2,:,:),4)))
xlim([0, 100])
title('100 & 300ms loss full view')
subplot(9,1,8)
plot(squeeze(mean(losses_100_300ms(:,2,:,:),4)))
ylim([0, 40000])
xlim([0, 100])
title('100 & 300ms zoomed')
subplot(9,1,9)
plot(squeeze(mean(mean(losses_100_300ms(:,2,:,:),4),3)))
ylim([0, 40000])
xlim([0, 100])
title('100 & 300ms mean_all')

% Get the best from each run
Rasters_100ms = zeros(220,10,29801);
PSTHs_100ms = zeros([220,149]);
[val_100,idx_100] = min(squeeze(losses_100ms(100,2,:,:))');
Rasters_200ms = zeros(220,10,29801);
PSTHs_200ms = zeros([220,149]);
[val_200,idx_200] = min(squeeze(losses_200ms(30,2,:,:))');
Rasters_100_300ms = zeros(220,10,29801);
PSTHs_100_300ms = zeros([220,149]);
[val_100_300,idx_100_300] = min(squeeze(losses_100_300ms(50,2,:,:))');


Rasters_data = zeros(220,10,29801);
PSTHs_data = zeros([220,149]);


bin_width = 200; %In 0.1 ms
bin_edges = 1:bin_width:29801+1;
dummy_idxs = 1:29801;

for k = 1:220
    Rasters_100ms(k,:,:) = squeeze(output_100ms(k,idx_100(k),:,:,:));
    PSTHs_100ms(k,:) = histcounts(squeeze(Rasters_100ms(k,:,:)).*dummy_idxs,bin_edges);
    Rasters_200ms(k,:,:) = squeeze(output_200ms(k,idx_200(k),:,:,:));
    PSTHs_200ms(k,:) = histcounts(squeeze(Rasters_200ms(k,:,:)).*dummy_idxs,bin_edges);
    Rasters_100_300ms(k,:,:) = squeeze(output_100_300ms(k,idx_100_300(k),:,:,:));
    PSTHs_100_300ms(k,:) = histcounts(squeeze(Rasters_100_300ms(k,:,:)).*dummy_idxs,bin_edges);
    
    SpikeTimes = all_data(k).ctrl_tar1_timestamps(:,1);
    picture = zeros(10,29801);
    for m = 1:10
        stim_mask = logical((SpikeTimes{m} > 0) .* (SpikeTimes{m} < 2.9801));
        trial_indicies = round(SpikeTimes{m}(stim_mask)*10000);
        picture(m,trial_indicies+1) = 1;
    end
    Rasters_data(k,:,:) = picture;
    PSTHs_data(k,:) = histcounts(squeeze(Rasters_data(k,:,:)).*dummy_idxs,bin_edges);
end

%%

%Loss comparison at 30epochs (Just L2 direct loss not PSTH)
x = 0;
y = 0;
figure;
subplot(1,2,1)
for k = 1:220
    plot(1,squeeze(losses_100ms(30,1,k,idx_100(k))),'k.'); hold on
    plot(2,squeeze(losses_200ms(30,1,k,idx_200(k))),'k.'); hold on
    if (losses_100ms(30,1,k,idx_100(k))) > squeeze(losses_200ms(30,1,k,idx_200(k)))
        plot([1,2],[squeeze(losses_100ms(30,1,k,idx_100(k))),squeeze(losses_200ms(30,1,k,idx_200(k)))],'b--'); hold on
        x = x+1;
    else
        plot([1,2],[squeeze(losses_100ms(30,1,k,idx_100(k))),squeeze(losses_200ms(30,1,k,idx_200(k)))],'r--'); hold on
        y = y+1;
    end
    
end
xticks([1 2])
xlim([0.25 2.75])
xticklabels({'100 ms','200 ms'})
title('loss comparison (L2 loss): ' + string(y) + ' out of 220 are better with 100ms loss')

x = 0;
y = 0;
subplot(1,2,2)
for k = 1:220
    plot(1,squeeze(losses_100ms(30,2,k,idx_100(k))),'k.'); hold on
    plot(2,squeeze(losses_200ms(30,2,k,idx_200(k))),'k.'); hold on
    if (losses_100ms(30,2,k,idx_100(k))) > squeeze(losses_200ms(30,2,k,idx_200(k)))
        plot([1,2],[squeeze(losses_100ms(30,2,k,idx_100(k))),squeeze(losses_200ms(30,2,k,idx_200(k)))],'b--'); hold on
        x = x+1;
    else
        plot([1,2],[squeeze(losses_100ms(30,2,k,idx_100(k))),squeeze(losses_200ms(30,2,k,idx_200(k)))],'r--'); hold on
        y = y+1;
    end
    
end
xticks([1 2])
xlim([0.25 2.75])
xticklabels({'100 ms','200 ms'})
title('loss comparison (200ms PSTH): ' + string(y) + ' out of 220 are better with 100ms loss')


% Lest look at the same units from before
%Cell 102, 133, 165, 200, 210, 84, 7, 36, 34, 79
cells = [102,133,165,200,210,84,7,36,34,79];

% Lets look at some of the worst fits in the 200ms case and see how the
% 100ms did

%n = numel(idx_200);
%out = arrayfun(@(k) losses_200ms(30,2,k,idx_200(k)), 1:n);

%[vals_200_worst, idx_200_worst] = maxk(out, 10);
%cells = idx_200_worst;

%Lets look at the worst by correlation
cors = [];
for k = 1:220
    cors = [cors, max(xcorr(PSTHs_100ms(k,:),PSTHs_data(k,:),'coeff'))];
end
[vals_200_worst, idx_200_worst] = maxk(cors, 10);
%cells = idx_200_worst;

figure;
for k = 1:10
    subplot(10,4,(k-1)*4+1)
    spy(squeeze(Rasters_100ms(cells(k),:,:)),'k')
    title('Cell ' + string(cells(k)) + ' 100ms model')
    subplot(10,4,(k-1)*4+2)
    spy(squeeze(Rasters_200ms(cells(k),:,:)),'k')
    title('Cell ' + string(cells(k)) + ' 200ms model')
    subplot(10,4,(k-1)*4+3)
    spy(squeeze(Rasters_data(cells(k),:,:)),'k')
    title('Cell ' + string(cells(k)) + ' Data')
    subplot(10,4,(k-1)*4+4)
    plot(squeeze(PSTHs_100ms(cells(k),:)),'b'); hold on
    plot(squeeze(PSTHs_200ms(cells(k),:)),'r'); hold on
    plot(squeeze(PSTHs_data(cells(k),:)),'k'); hold on
    %title('Cell ' + string(cells(k)) + ' PSTH comparison. Xorr: ' + string(vals_200_worst(k)))
end

%%
cors_100 = [];
cors_200 = [];
cors_100_300 = [];
for k = 1:220
    cors_100 = [cors_100, max(xcorr(PSTHs_100ms(k,:),PSTHs_data(k,:),'coeff'))];
    cors_200 = [cors_200, max(xcorr(PSTHs_200ms(k,:),PSTHs_data(k,:),'coeff'))];
    cors_100_300 = [cors_100_300, max(xcorr(PSTHs_100_300ms(k,:),PSTHs_data(k,:),'coeff'))];
end

%Lets try ploting every single cell
cells_per_page = 20;
n_cells = 220;
n_pages = ceil(n_cells / cells_per_page);

for p = 1:n_pages
    start_k = (p-1)*cells_per_page + 1;
    end_k   = min(p*cells_per_page, n_cells);
    n_rows_this_page = end_k - start_k + 1;

    f = figure('Visible','off', ...
               'Units','pixels', ...
               'Position',[100 100 2600 5000], ...
               'Color','w');

    t = tiledlayout(n_rows_this_page,5, ...
                    'TileSpacing','tight', ...
                    'Padding','tight');

    for k = start_k:end_k

        nexttile
        A = squeeze(Rasters_100ms(k,:,:));
        plotSpikeRaster(A)
        title("Cell " + k + " 100ms model "  + "Xcorr: " + num2str(cors_100(k),'%.2f'),'FontSize',7)
        
        nexttile
        A = squeeze(Rasters_200ms(k,:,:));
        plotSpikeRaster(A)
        title("Cell " + k + " 200ms model " + "Xcorr: " + num2str(cors_200(k),'%.2f'),'FontSize',7)
        
        nexttile
        A = squeeze(Rasters_100_300ms(k,:,:));
        plotSpikeRaster(A)
        title("Cell " + k + " 100 & 300ms model "+ "Xcorr: " + num2str(cors_100_300(k),'%.2f'),'FontSize',7)
        
        nexttile
        A = squeeze(Rasters_data(k,:,:));
        plotSpikeRaster(A)
        title("Cell " + k + " Data",'FontSize',7)

        % --- PSTH ---
        nexttile
        plot(squeeze(PSTHs_100ms(k,:)),'b','LineWidth',0.75); hold on
        plot(squeeze(PSTHs_200ms(k,:)),'r','LineWidth',0.75)
        plot(squeeze(PSTHs_100_300ms(k,:)),'g','LineWidth',0.75)
        plot(squeeze(PSTHs_data(k,:)),'k','LineWidth',0.75)
        hold off
        set(gca,'XTick',[],'Box','on')
        title("Cell " + k + " PSTH. ",'FontSize',7)
    end

    % optional: remove this too if you want a little more room
    sgtitle("Cells " + start_k + " to " + end_k, ...
            'FontSize',12,'FontWeight','bold')

    exportgraphics(f, sprintf('all_cells_page_%02d.pdf', p), 'ContentType','vector')
    exportgraphics(f, sprintf('all_cells_page_%02d.png', p), 'Resolution',300)

    close(f)
end



%% Look at noise corrected correlation distributions between loss values


bin_width = 200; %In 0.1 ms

%Rasters_100ms_trial_averaged = zeros(220,10,29801);
PSTHs_data_trial_averaged = zeros([220,floor(29801/bin_width)]);


bin_edges = 1:bin_width:29801+1;
dummy_idxs = 1:29801;

for k = 1:220
    %Rasters_100ms_trial_averaged(k,:,:) = squeeze(output_100ms(k,idx_100(k),:,:,:));
    PSTHs_data_trial_averaged(k,:) = histcounts(squeeze(Rasters_data(k,:,:)).*dummy_idxs,bin_edges)/10; %10 trials
end

PSTHs_data_varience = zeros([220,floor(29801/bin_width)]);
PSTHs_data_total_varience = zeros([220,1]);
PSTHs_data_varience_out = zeros([220,1]);

for z = 1:220
    PSTHs_data_total_varience(z) = (1/((floor(29801/bin_width))-1))*sum((PSTHs_data_trial_averaged(z,:) - (1/floor(29801/bin_width))*sum(PSTHs_data_trial_averaged(z,:))).^2);
    for k = 1:floor(29801/bin_width)
        for m = 1:10
            PSTHs_data_varience(z,k) =  PSTHs_data_varience(z,k) + (sum(Rasters_data(z,m,(1+(k-1)*220):(k)*220)) - PSTHs_data_trial_averaged(z,k)).^2;
        end
        PSTHs_data_varience(z,k) = PSTHs_data_varience(z,k)/9;
    end
    PSTHs_data_varience_out(z) = sum(PSTHs_data_varience(z,:))/floor(29801/bin_width);
end


signal_noise = PSTHs_data_total_varience - PSTHs_data_varience_out/10;

covs = zeros([220,1]);
pncs = zeros([220,1]);
for k = 1:220
    mean_sub_PSTH = histcounts(squeeze(Rasters_100ms(k,:,:)).*dummy_idxs,bin_edges) - (1/floor(29801/bin_width))*sum(histcounts(squeeze(Rasters_100ms(k,:,:)).*dummy_idxs,bin_edges));
    mean_sub_PSTH_data = histcounts(squeeze(Rasters_data(k,:,:)).*dummy_idxs,bin_edges) - (1/floor(29801/bin_width))*sum(histcounts(squeeze(Rasters_data(k,:,:)).*dummy_idxs,bin_edges));
    sigmap = sum(mean_sub_PSTH.^2)/(floor(29801/bin_width)-1);
    covs(k)  = sum(mean_sub_PSTH .* mean_sub_PSTH_data) / (numel(mean_sub_PSTH) - 1);
    pncs(k) = covs(k)/sqrt(signal_noise(k)*sigmap);
end





%%
bin_width = 100; % in 0.1 ms samples

[n_neurons, n_trials, n_time] = size(Rasters_data);
n_bins = floor(n_time / bin_width);   % drops leftover samples at the end, if any
n_time_used = n_bins * bin_width;

PSTHs_data_trial_averaged = zeros(n_neurons, n_bins);
PSTHs_model_trial_averaged = zeros(n_neurons, n_bins);

PSTHs_data_varience = zeros(n_neurons, n_bins);
PSTHs_data_total_varience = zeros(n_neurons, 1);
PSTHs_data_varience_out = zeros(n_neurons, 1);

signal_noise = zeros(n_neurons, 1);
covs = zeros(n_neurons, 1);
pncs = nan(n_neurons, 1);

for z = 1:n_neurons
    data_bin_counts = zeros(n_trials, n_bins);
    model_bin_counts = zeros(n_trials, n_bins);

    for m = 1:n_trials
        data_trial = squeeze(Rasters_data(z, m, 1:n_time_used));
        model_trial = squeeze(Rasters_100ms(z, m, 1:n_time_used));

        data_trial = reshape(data_trial, bin_width, n_bins);
        model_trial = reshape(model_trial, bin_width, n_bins);

        data_bin_counts(m, :) = sum(data_trial, 1);
        model_bin_counts(m, :) = sum(model_trial, 1);
    end

    % Trial-averaged PSTHs
    PSTHs_data_trial_averaged(z, :) = mean(data_bin_counts, 1);
    PSTHs_model_trial_averaged(z, :) = mean(model_bin_counts, 1);

    % Total variance of the data PSTH across time bins
    PSTHs_data_total_varience(z) = var(PSTHs_data_trial_averaged(z, :), 0, 2);

    % Trial-to-trial noise variance, averaged across bins
    PSTHs_data_varience(z, :) = var(data_bin_counts, 0, 1);
    PSTHs_data_varience_out(z) = mean(PSTHs_data_varience(z, :));

    % Signal variance estimate
    signal_noise(z) = PSTHs_data_total_varience(z) - PSTHs_data_varience_out(z) / n_trials;

    % Mean-subtracted PSTHs
    mean_sub_PSTH = PSTHs_model_trial_averaged(z, :) - mean(PSTHs_model_trial_averaged(z, :));
    mean_sub_PSTH_data = PSTHs_data_trial_averaged(z, :) - mean(PSTHs_data_trial_averaged(z, :));

    % Model variance and data-model covariance
    sigmap = sum(mean_sub_PSTH .^ 2) / (n_bins - 1);
    covs(z) = sum(mean_sub_PSTH .* mean_sub_PSTH_data) / (n_bins - 1);

    % Noise-corrected correlation
    if signal_noise(z) > 0 && sigmap > 0
        pncs(z) = covs(z) / sqrt(signal_noise(z) * sigmap);
    end
end




figure;
histogram(pncs)

figure;
histogram(cors_100)

%%
bin_width = 100; %In 0.1 ms

Rasters_data = zeros(220,10,29801);
PSTHs_data_100 = zeros([220,floor(29801/bin_width)]);

bin_edges = 1:bin_width:29801+1;
dummy_idxs = 1:29801;

PSTHs_100ms = zeros([220,298]);
for k = 1:220
    PSTHs_100ms(k,:) = histcounts(squeeze(Rasters_100ms(k,:,:)).*dummy_idxs,bin_edges);
    SpikeTimes = all_data(k).ctrl_tar1_timestamps(:,1);
    picture = zeros(10,29801);
    for m = 1:10
        stim_mask = logical((SpikeTimes{m} > 0) .* (SpikeTimes{m} < 2.9801));
        trial_indicies = round(SpikeTimes{m}(stim_mask)*10000);
        picture(m,trial_indicies+1) = 1;
    end
    Rasters_data(k,:,:) = picture;
    PSTHs_data_100(k,:) = histcounts(squeeze(Rasters_data(k,:,:)).*dummy_idxs,bin_edges);
end

bin_width = 200; %In 0.1 ms

Rasters_data = zeros(220,10,29801);
PSTHs_data_200 = zeros([220,floor(29801/bin_width)]);

bin_edges = 1:bin_width:29801+1;
dummy_idxs = 1:29801;



for k = 1:220
    
    SpikeTimes = all_data(k).ctrl_tar1_timestamps(:,1);
    picture = zeros(10,29801);
    for m = 1:10
        stim_mask = logical((SpikeTimes{m} > 0) .* (SpikeTimes{m} < 2.9801));
        trial_indicies = round(SpikeTimes{m}(stim_mask)*10000);
        picture(m,trial_indicies+1) = 1;
    end
    Rasters_data(k,:,:) = picture;
    PSTHs_data_200(k,:) = histcounts(squeeze(Rasters_data(k,:,:)).*dummy_idxs,bin_edges);
end



cors_100 = [];
cors_200 = [];
cors_100_300 = [];
cors_100_300_200 = [];
for k = 1:220
    cors_100 = [cors_100, max(xcorr(PSTHs_100ms(k,:),PSTHs_data_100(k,:),'coeff'))];
    cors_200 = [cors_200, max(xcorr(PSTHs_200ms(k,:),PSTHs_data_200(k,:),'coeff'))];
    %cors_100_300 = [cors_100_300, max(xcorr(PSTHs_100_300ms(k,:),PSTHs_data_100(k,:),'coeff'))];
    cors_100_300_200 = [cors_100_300_200, max(xcorr(PSTHs_100_300ms(k,:),PSTHs_data_200(k,:),'coeff'))];
end

figure(Position=[0,0,1000,1000]);
subplot(2,2,1);
histogram(cors_100)
subplot(2,2,2);
histogram(cors_200)
subplot(2,2,3);
%histogram(cors_100_300)
subplot(2,2,4);
histogram(cors_100_300_200)





%% Statistical different calculation

hs = [];
ps = [];
for k = 1:220
    
    pre_vals = [];

    for m = 1:10
        for z  = 1:4
            SpikeTimes = all_data(k).ctrl_tar1_timestamps(m,z);
            pre_vals = [pre_vals, length(SpikeTimes{1})];
        end
    end


    post_vals = [];

    for m = 1:10
        for z  = 1:4
            SpikeTimes = all_data(k).laser_tar1_timestamps(m,z);
            post_vals = [post_vals, length(SpikeTimes{1})];
        end
    end

    %Seems appropriate to do a paired t-test since each fr has a
    %corresponding fr post-laser

    %H0 is that post_vals < pre_vals
    %diffs = pre_vals-post_vals;
    %t = mean(diffs)/(std(diffs)/sqrt(40));
    [h,p] = ttest(pre_vals,post_vals,"Tail","right");
    

    hs = [hs,h];
    ps = [ps,p];

end

figure;
stem(ps); hold on
ylim([0, 0.1])






%% Replot the cells but with PV and SST cells excluded and the new run and the data iwth the peaks of the tuning curves.
load('all_units_info_with_polished_criteria_modified_perf.mat','all_data');
load('output_compressed_Eprop_All_cells_20260402_090929')

bin_width = 100; %In 0.1 ms

Rasters_data = zeros(220,10,29801);
PSTHs_data= zeros([220,floor(29801/bin_width)]);

bin_edges = 1:bin_width:29801+1;
dummy_idxs = 1:29801;

PSTHs_sim = zeros([220,floor(29801/bin_width)]);
Rasters_sim = zeros(220,10,29801);

[val,idx] = min(squeeze(losses(50,2,:,:))');

Nss = [];

for k = 1:220
    %Test if it is narrow spiking
    Nss = [Nss, all_data(k).is_NS];

    if  all_data(k).is_NS == 0

        %if hs(k) == 0

            %Test if there is significant supression by the laser
            %Look at the peak response pre-laser, look at it post laser,
            %compare
    
    
            Rasters_sim(k,:,:) = squeeze(output(k,idx(k),:,:,:));
            PSTHs_sim(k,:) = histcounts(squeeze(Rasters_sim(k,:,:)).*dummy_idxs,bin_edges);
            SpikeTimes = all_data(k).ctrl_tar1_timestamps(:,1);
            picture = zeros(10,29801);
            for m = 1:10
                stim_mask = logical((SpikeTimes{m} > 0) .* (SpikeTimes{m} < 2.9801));
                trial_indicies = round(SpikeTimes{m}(stim_mask)*10000);
                picture(m,trial_indicies+1) = 1;
            end
            Rasters_data(k,:,:) = picture;
            PSTHs_data(k,:) = histcounts(squeeze(Rasters_data(k,:,:)).*dummy_idxs,bin_edges);
        %end
    end
end


mask = Nss;
mask(mask>1) = 1;





%% Direct FR comparison

pre_laser = [];
post_laser = [];
for k = 1:220

    tuning = all_data(k).tuning_type;
    focus = 0;

    if strcmp(tuning, 'contra-tuned')
        focus = 1;
    elseif strcmp(tuning, '45°-tuned')
        focus = 2;
    elseif strcmp(tuning, 'center-tuned')
        focus = 3;
    elseif strcmp(tuning, 'ipsi-tuned')
        focus = 4;
    end

    spk_val_pre = 0;
    SpikeTimes = all_data(k).ctrl_tar1_timestamps(:,focus);
    for m = 1:10
        spk_val_pre = spk_val_pre + length(SpikeTimes{m});
    end
    pre_laser = [pre_laser, spk_val_pre];

    spk_val_post = 0;
    SpikeTimes = all_data(k).laser_tar1_timestamps(:,focus);
    for m = 1:10
        spk_val_post = spk_val_post + length(SpikeTimes{m});
    end
    post_laser = [post_laser, spk_val_post];



end

perc_decrease = (pre_laser-post_laser)./pre_laser;

figure;
n = 1:220;
plot(n, perc_decrease)

figure;
n = 1:220;
plot(n, perc_decrease)
ylim([0,1])

%%

%Lets try ploting every single cell
cells_per_page = 20;
n_cells = 220-sum(mask);
n_pages = ceil(n_cells / cells_per_page);

%Mask everything
cells = 1:220;
cells = cells(mask==0);

Rasters_sim = Rasters_sim(mask==0,:,:);
Rasters_data = Rasters_data(mask==0,:,:);

PSTHs_sim = PSTHs_sim(mask==0,:,:);
PSTHs_data = PSTHs_data(mask==0,:,:);

for p = 1:n_pages
    start_k = (p-1)*cells_per_page + 1;
    end_k   = min(p*cells_per_page, n_cells);
    n_rows_this_page = end_k - start_k + 1;

    f = figure('Visible','off', ...
               'Units','pixels', ...
               'Position',[100 100 2600 5000], ...
               'Color','w');

    t = tiledlayout(n_rows_this_page,3, ...
                    'TileSpacing','tight', ...
                    'Padding','tight');

    for k = start_k:end_k

        nexttile
        A = squeeze(Rasters_sim(k,:,:));
        plotSpikeRaster(A)
        title("Cell " + cells(k) + " 10/30ms model "  + "Xcorr: ",'FontSize',7)
        
        nexttile
        A = squeeze(Rasters_data(k,:,:));
        plotSpikeRaster(A)
        title("Cell " + cells(k) + " data " + "Xcorr: ",'FontSize',7)

        % --- PSTH ---
        nexttile
        plot(squeeze(PSTHs_sim(k,:)),'g','LineWidth',0.75); hold on
        plot(squeeze(PSTHs_data(k,:)),'k','LineWidth',0.75)
        hold off
        set(gca,'XTick',[],'Box','on')
        title("Cell " + cells(k) + " PSTH. ",'FontSize',7)
    end

    % optional: remove this too if you want a little more room
    sgtitle("Cells " + cells(start_k) + " to " + cells(end_k), ...
            'FontSize',12,'FontWeight','bold')

    exportgraphics(f, sprintf('all_cells_page_%02d.pdf', p), 'ContentType','vector')
    exportgraphics(f, sprintf('all_cells_page_%02d.png', p), 'Resolution',300)

    close(f)
end




%% MDS plots

num_cells2 = 198;

comp_matrix_data = zeros([num_cells2,num_cells2]);
comp_matrix_sim = zeros([num_cells2,num_cells2]);
for m = 1:num_cells2
    for k = 1:num_cells2
        comp_matrix_data(m,k) = sqrt(sum((PSTHs_data(m,:) - PSTHs_data(k,:)).^2));  %RMSE
        comp_matrix_sim(m,k) = sqrt(sum((PSTHs_sim(m,:) - PSTHs_sim(k,:)).^2));  
    end
end

[Y_data,stress_data] = mdscale(comp_matrix_data,2, 'criterion','metricsstress');
[Y_sim,stress_sim] = mdscale(comp_matrix_sim,2, 'criterion','metricsstress');

% Procruste
[d,Z,transform] = procrustes(Y_data,Y_sim);


%Transformation and procruste without log transform

figure;
scatter(Y_data(:,1),Y_data(:,2),36,'ro','filled'); hold on
labels = string(cells);
text(Y_data(:,1),Y_data(:,2),labels,'FontSize',10)
%set(gca,'XScale','log','YScale','log') %Negetive values. Can still do log
%compression but will need to finagle it.
xlabel('mds axis 1')
ylabel('mds axis 2')

figure;
scatter(Y_sim(:,1),Y_sim(:,2),36,'bo','filled'); hold on
labels = string(cells);
%text(Y_sim(:,1),Y_sim(:,2),labels,'FontSize',10)
%set(gca,'XScale','log','YScale','log') %Negetive values. Can still do log
%compression but will need to finagle it.
xlabel('mds axis 1')
ylabel('mds axis 2')

figure;
scatter(Y_data(:,1),Y_data(:,2),36,'ro','filled'); hold on
labels = string(cells);
text(Y_data(:,1),Y_data(:,2),labels,'FontSize',10,'Color',[1,0,0])

scatter(Z(:,1),Z(:,2),36,'ko','filled');
text(Z(:,1),Z(:,2),labels,'FontSize',10)

dists = [];

for k = 1:184
    plot([Y_data(k,1),Z(k,1)],[Y_data(k,2),Z(k,2)], 'Color', [1 0 0 0.1])
    dists = [dists, sqrt((Y_data(k,1)-Z(k,1))^2 + (Y_data(k,2)-Z(k,2))^2)];
end

%set(gca,'XScale','log','YScale','log') %Negetive values. Can still do log
%compression but will need to finagle it.
xlabel('mds axis 1')
ylabel('mds axis 2')





layer = [];

for k = 1:220
    if  strcmp(all_data(k).layer,'L2/3')
        layer = [layer,1];
    elseif strcmp(all_data(k).layer, 'L5/6')
        layer = [layer,2];
    elseif strcmp(all_data(k).layer,'L4')
        layer = [layer,3];
    elseif strcmp(all_data(k).layer, 'NaN')
        layer = [layer,4];
    end
end

layer = layer(mask==0);

figure(position=[0,0,1000,1000]);%Data plot

%Push all data beyond x & y == 1
Y_data(:,1) = (Y_data(:,1) - min(Y_data(:,1)))+1;
Y_data(:,2) = (Y_data(:,2) - min(Y_data(:,2)))+1;

xplot = log(Y_data(:,1));
yplot = log(Y_data(:,2));

scatter(xplot,yplot,36,layer,'filled'); hold on
labels = string(cells);
text(xplot,yplot,labels,'FontSize',10)
%set(gca,'XScale','log','YScale','log') %Negetive values. Can still do log
%compression but will need to finagle it.
xlabel('mds axis 1')
ylabel('mds axis 2')

figure(position=[0,0,1000,1000]); %Sim plot
Y_sim(:,1) = (Y_sim(:,1) - min(Y_sim(:,1)))+1;
Y_sim(:,2) = (Y_sim(:,2) - min(Y_sim(:,2)))+1;

xplot = log(Y_sim(:,1));
yplot = log(Y_sim(:,2));
scatter(xplot,yplot,36,layer,'filled'); hold on
labels = string(cells);
text(xplot,yplot,labels,'FontSize',10)
xlabel('mds axis 1')
ylabel('mds axis 2')

figure(position=[0,0,1000,1000]); %Sim plot (y-axis flipped)
Y_sim(:,1) = (Y_sim(:,1) - min(Y_sim(:,1)))+1;
Y_sim(:,2) = (Y_sim(:,2) - min(Y_sim(:,2)))+1;

xplot = log(Y_sim(:,1));
yplot = log(Y_sim(:,2));
scatter(xplot,-yplot,36,layer,'filled'); hold on
labels = string(cells);
text(xplot,-yplot,labels,'FontSize',10)
xlabel('mds axis 1')
ylabel('mds axis 2')


% Plot along the branches


%% MDS Axis 1 plot


f = figure('Visible','on', ...
               'Units','pixels', ...
               'Position',[100 100 2600 5000], ...
               'Color','w');

t = tiledlayout(10,3, ...
                'TileSpacing','tight', ...
                'Padding','tight');

%cells_along_mds_axis1 = [206,160,166,208,69,133,193,111,42,1];
%Note - this is just looking at some of the poor fits according to the
%noise corrected correlation
noise_corrected_bad = [37,156,55,185,51,62,17,113,69,173];
cells_along_mds_axis1 = noise_corrected_bad;

cells_back = [];

for k = 1:10
    
    cells_back = [cells_back, find(cells == cells_along_mds_axis1(k))];
    
end

for k = 1:10

    nexttile
    A = squeeze(Rasters_sim(cells_back(k),:,:));
    plotSpikeRaster(A)
    title("Cell " + cells_back(k) + " 10/30ms model "  + "Xcorr: ",'FontSize',7)
    
    nexttile
    A = squeeze(Rasters_data(cells_back(k),:,:));
    plotSpikeRaster(A)
    title("Cell " + cells_back(k) + " data " + "Xcorr: ",'FontSize',7)

    % --- PSTH ---
    nexttile
    plot(squeeze(PSTHs_sim(cells_back(k),:)),'g','LineWidth',0.75); hold on
    plot(squeeze(PSTHs_data(cells_back(k),:)),'k','LineWidth',0.75)
    hold off
    set(gca,'XTick',[],'Box','on')
    title("Cell " + cells_back(k) + " PSTH. ",'FontSize',7)
end

% optional: remove this too if you want a little more room
sgtitle("MDS Axis 1")

exportgraphics(f, sprintf('MDS_axis1%02d.pdf', p), 'ContentType','vector')

%% Color MDS by rate

data_rate = [];

for k = 1:num_cells2
    data_rate = [data_rate,sum(sum(Rasters_data(k,:,:)))/2.98/10];
end

model_rate = [];

for k = 1:num_cells2
    model_rate = [model_rate,sum(sum(Rasters_sim(k,:,:)))/2.98/10];
end


figure;
Y_data(:,1) = (Y_data(:,1) - min(Y_data(:,1)))+1;
Y_data(:,2) = (Y_data(:,2) - min(Y_data(:,2)))+1;

xplot = log(Y_data(:,1));
yplot = log(Y_data(:,2));
scatter(xplot,yplot,36,data_rate,'filled'); hold on
labels = string(cells);
text(xplot,yplot,labels,'FontSize',10)
%set(gca,'XScale','log','YScale','log') %Negetive values. Can still do log
%compression but will need to finagle it.
xlabel('mds axis 1')
ylabel('mds axis 2')

figure(position=[0,0,1000,1000]); %Sim plot (y-axis flipped)
Y_sim(:,1) = (Y_sim(:,1) - min(Y_sim(:,1)))+1;
Y_sim(:,2) = (Y_sim(:,2) - min(Y_sim(:,2)))+1;

xplot = log(Y_sim(:,1));
yplot = log(Y_sim(:,2));
scatter(xplot,-yplot,36,model_rate,'filled'); hold on
labels = string(cells);
text(xplot,-yplot,labels,'FontSize',10)
xlabel('mds axis 1')
ylabel('mds axis 2')


%% MDS Axis 2 plot


f = figure('Visible','on', ...
               'Units','pixels', ...
               'Position',[100 100 2600 5000], ...
               'Color','w');

t = tiledlayout(10,3, ...
                'TileSpacing','tight', ...
                'Padding','tight');

cells_along_mds_axis2 = [106,6,78,27,53,202,210,186,200,7];

cells_back = [];

for k = 1:10
    
    cells_back = [cells_back, find(cells == cells_along_mds_axis2(k))];
    
end

for k = 1:10

    nexttile
    A = squeeze(Rasters_sim(cells_back(k),:,:));
    plotSpikeRaster(A)
    title("Cell " + cells_back(k) + " 10/30ms model "  + "Xcorr: ",'FontSize',7)
    
    nexttile
    A = squeeze(Rasters_data(cells_back(k),:,:));
    plotSpikeRaster(A)
    title("Cell " + cells_back(k) + " data " + "Xcorr: ",'FontSize',7)

    % --- PSTH ---
    nexttile
    plot(squeeze(PSTHs_sim(cells_back(k),:)),'g','LineWidth',0.75); hold on
    plot(squeeze(PSTHs_data(cells_back(k),:)),'k','LineWidth',0.75)
    hold off
    set(gca,'XTick',[],'Box','on')
    title("Cell " + cells_back(k) + " PSTH. ",'FontSize',7)
end

% optional: remove this too if you want a little more room
sgtitle("MDS Axis 2")

exportgraphics(f, sprintf('MDS_axis2%02d.pdf', p), 'ContentType','vector')


%% Three Branches


f = figure('Visible','on', ...
               'Units','pixels', ...
               'Position',[100 100 2600 5000], ...
               'Color','w');

t = tiledlayout(10,3, ...
                'TileSpacing','tight', ...
                'Padding','tight');

cells_along_mds_b1 = [206,24,66,151,105,161,166,80,170,40];

cells_back = [];

for k = 1:10
    
    cells_back = [cells_back, find(cells == cells_along_mds_b1(k))];
    
end

for k = 1:10

    nexttile
    A = squeeze(Rasters_sim(cells_back(k),:,:));
    plotSpikeRaster(A)
    title("Cell " + cells_back(k) + " 10/30ms model "  + "Xcorr: ",'FontSize',7)
    
    nexttile
    A = squeeze(Rasters_data(cells_back(k),:,:));
    plotSpikeRaster(A)
    title("Cell " + cells_back(k) + " data " + "Xcorr: ",'FontSize',7)

    % --- PSTH ---
    nexttile
    plot(squeeze(PSTHs_sim(cells_back(k),:)),'g','LineWidth',0.75); hold on
    plot(squeeze(PSTHs_data(cells_back(k),:)),'k','LineWidth',0.75)
    hold off
    set(gca,'XTick',[],'Box','on')
    title("Cell " + cells_back(k) + " PSTH. ",'FontSize',7)
end

% optional: remove this too if you want a little more room
sgtitle("Branch 1")

exportgraphics(f, sprintf('Branch 1%02d.pdf', p), 'ContentType','vector')


%% Branch 2 (Upper)

f = figure('Visible','on', ...
               'Units','pixels', ...
               'Position',[100 100 2600 5000], ...
               'Color','w');

t = tiledlayout(10,3, ...
                'TileSpacing','tight', ...
                'Padding','tight');

cells_along_mds_b1 = [133,102,106,9,12,122,6,12,28,157];

cells_back = [];

for k = 1:10
    
    cells_back = [cells_back, find(cells == cells_along_mds_b1(k))];
    
end

for k = 1:10

    nexttile
    A = squeeze(Rasters_sim(cells_back(k),:,:));
    plotSpikeRaster(A)
    title("Cell " + cells_back(k) + " 10/30ms model "  + "Xcorr: ",'FontSize',7)
    
    nexttile
    A = squeeze(Rasters_data(cells_back(k),:,:));
    plotSpikeRaster(A)
    title("Cell " + cells_back(k) + " data " + "Xcorr: ",'FontSize',7)

    % --- PSTH ---
    nexttile
    plot(squeeze(PSTHs_sim(cells_back(k),:)),'g','LineWidth',0.75); hold on
    plot(squeeze(PSTHs_data(cells_back(k),:)),'k','LineWidth',0.75)
    hold off
    set(gca,'XTick',[],'Box','on')
    title("Cell " + cells_back(k) + " PSTH. ",'FontSize',7)
end

% optional: remove this too if you want a little more room
sgtitle("Branch 2")

exportgraphics(f, sprintf('Branch 2%02d.pdf', p), 'ContentType','vector')

%% Branch 3 (Lower)

f = figure('Visible','on', ...
               'Units','pixels', ...
               'Position',[100 100 2600 5000], ...
               'Color','w');

t = tiledlayout(10,3, ...
                'TileSpacing','tight', ...
                'Padding','tight');

cells_along_mds_b1 = [7,200,210,68,191,132,13,202,107,124];

cells_back = [];

for k = 1:10
    
    cells_back = [cells_back, find(cells == cells_along_mds_b1(k))];
    
end

for k = 1:10

    nexttile
    A = squeeze(Rasters_sim(cells_back(k),:,:));
    plotSpikeRaster(A)
    title("Cell " + cells_back(k) + " 10/30ms model "  + "Xcorr: ",'FontSize',7)
    
    nexttile
    A = squeeze(Rasters_data(cells_back(k),:,:));
    plotSpikeRaster(A)
    title("Cell " + cells_back(k) + " data " + "Xcorr: ",'FontSize',7)

    % --- PSTH ---
    nexttile
    plot(squeeze(PSTHs_sim(cells_back(k),:)),'g','LineWidth',0.75); hold on
    plot(squeeze(PSTHs_data(cells_back(k),:)),'k','LineWidth',0.75)
    hold off
    set(gca,'XTick',[],'Box','on')
    title("Cell " + cells_back(k) + " PSTH. ",'FontSize',7)
end

% optional: remove this too if you want a little more room
sgtitle("Branch 3")

exportgraphics(f, sprintf('Branch 3%02d.pdf', p), 'ContentType','vector')

%% Maybe axis 2 is first ttp latency?

ttp_data = [];
ttp_model = [];

%Quantify
for k = 1:num_cells2
    
    %Look in the first 320ms
    [val,idx] = max(PSTHs_data(k,1:16));
    ttp_data = [ttp_data, idx];
    [val,idx] = max(PSTHs_sim(k,1:16));
    ttp_model = [ttp_model, idx];

end

figure;
Y_data(:,1) = (Y_data(:,1) - min(Y_data(:,1)))+1;
Y_data(:,2) = (Y_data(:,2) - min(Y_data(:,2)))+1;

xplot = log(Y_data(:,1));
yplot = log(Y_data(:,2));
scatter(xplot,yplot,36,ttp_data,'filled'); hold on
labels = string(cells);
text(xplot,yplot,labels,'FontSize',10)
%set(gca,'XScale','log','YScale','log') %Negetive values. Can still do log
%compression but will need to finagle it.
xlabel('mds axis 1')
ylabel('mds axis 2')

figure(position=[0,0,1000,1000]); %Sim plot (y-axis flipped)
Y_sim(:,1) = (Y_sim(:,1) - min(Y_sim(:,1)))+1;
Y_sim(:,2) = (Y_sim(:,2) - min(Y_sim(:,2)))+1;

xplot = log(Y_sim(:,1));
yplot = log(Y_sim(:,2));
scatter(xplot,-yplot,36,ttp_model,'filled'); hold on
labels = string(cells);
text(xplot,-yplot,labels,'FontSize',10)
xlabel('mds axis 1')
ylabel('mds axis 2')


function plotSpikeRaster(A)
    % A is trials x timebins, with 0/1 entries
    [trial, bin] = find(A > 0);

    if ~isempty(trial)
        x = [bin.'; bin.'; nan(1,numel(bin))];
        y = [trial.' - 0.4; trial.' + 0.4; nan(1,numel(trial))];
        line(x(:), y(:), 'Color', 'k', 'LineWidth', 0.5);
    end

    xlim([0.5, size(A,2)+0.5]);
    ylim([0.5, size(A,1)+0.5]);
    set(gca, ...
        'YDir','reverse', ...   % keeps trial 1 at the top like imagesc
        'XTick',[], ...
        'YTick',[], ...
        'Box','on');
end


