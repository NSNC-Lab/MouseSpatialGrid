%Loss over epochs visualization

%Avg loss
mean_over_batch = squeeze(mean(losses(:,2,:,:),4));
mean_over_all_cells = squeeze(mean(mean_over_batch,2));


figure;
subplot(2,1,1)
plot(mean_over_batch)
title('All cells losses')
ylabel('Loss')

subplot(2,1,2)
plot(mean_over_all_cells)
title('Average loss over all cells')
xlabel('epochs')
ylabel('Loss')

%%

%Compare 2 disparate FR cells (Cell 1 and 3)
%Note spontaneous FR was set wrong. Need to rerun and fix

%Get the params that we want to look at (last epoch, all params, cell 1 and
%3, all instances in batch)

all_data = [];

cell_data_1 = squeeze(params(end,:,1,:));
cell_data_3 = squeeze(params(end,:,117,:));

all_data = [cell_data_1,cell_data_3]';


[coeff, score, latent, tsquared, explained] = pca(all_data);


figure;
scatter(score(1:40,1),  score(1:40,2), 36, 'filled'); hold on
scatter(score(41:80,1), score(41:80,2), 36, 'filled'); hold on

scatter(mean(score(1:40,1)),  mean(score(1:40,2)), 120,"blue", 'filled'); hold on
scatter(mean(score(41:80,1)), mean(score(41:80,2)), 120,[1, 0.6, 0], 'filled');

dist = sqrt((mean(score(41:80,1))-mean(score(1:40,1)))^2 + (mean(score(41:80,2))-mean(score(1:40,2)))^2);


%%

%Find the cell that is "least like cell 1 by centroid distance"

all_dists = [];

for cells = 1:220
    all_data = [];
    
    cell_data_1 = squeeze(params(end,:,1,:));
    cell_data_3 = squeeze(params(end,:,cells,:));
    
    all_data = [cell_data_1,cell_data_3]';
    
    
    [coeff, score, latent, tsquared, explained] = pca(all_data);
    
    
    %figure;
    %scatter(score(1:40,1),  score(1:40,2), 36, 'filled'); hold on
    %scatter(score(41:80,1), score(41:80,2), 36, 'filled'); hold on
    
    %scatter(mean(score(1:40,1)),  mean(score(1:40,2)), 120,"blue", 'filled'); hold on
    %scatter(mean(score(41:80,1)), mean(score(41:80,2)), 120,[1, 0.6, 0], 'filled');
    
    dist = sqrt((mean(score(41:80,1))-mean(score(1:40,1)))^2 + (mean(score(41:80,2))-mean(score(1:40,2)))^2);
    all_dists = [all_dists,dist];
end



%%
%mu = mean(all_data(1,:), 2);
%sigma = std(all_data(1,:), 0, 2);
%Z = (all_data(1,:) - mu) ./ sigma;


%%
all_data = [];

for k = 1:220
    
    all_data = [all_data,squeeze(params(end,:,k,:))];

end

all_data_z_scored = zscore(all_data,0,2);
%cell_data_1 = squeeze(params(end,:,1,:));
%cell_data_3 = squeeze(params(end,:,117,:));

%all_data = reshape(squeeze(params(end,:,:,:)),[5,30,220]);
%all_data_reshaped = reshape(all_data,[5,220*30])';


[coeff, score, latent, tsquared, explained] = pca(all_data_z_scored');
figure('position',[0,0,500,500]);
scatter(score(:,1),  score(:,2), 36, 'filled'); hold on
%scatter(score(151:180,1),  score(151:180,2), 36, 'filled'); hold on
scatter(score(1:30,1),  score(1:30,2), 36, 'filled'); hold on
scatter(score(181:210,1),  score(181:210,2), 36, 'filled')

%%
for k = 1:220
    figure;
    scatter(score(30*(k-1)+1:30*(k),1),  score(30*(k-1)+1:30*(k),2), 36, 'filled');
end


%% Centroid Plot
reshape_score_data = reshape(score,[5,30,220]);
reshape_meaned_score_data = mean(reshape_score_data,2);

figure;
scatter(reshape_meaned_score_data(1,:),reshape_meaned_score_data(2,:))

%%



mu = [0 0];
Sigma = [1 0.6; 0.6 2];

mu2 = [10,10];

N = 5000;
S = mvnrnd(mu, Sigma, N);
S2 = mvnrnd(mu2, Sigma, N);

figure;
scatter(S(:,1), S(:,2), 6, '.'); axis equal; hold on;
scatter(S2(:,1), S2(:,2), 6, '.'); axis equal;
title('Samples from 2D Gaussian');


S_all = [S;S2];
[coeff, score, latent, tsquared, explained] = pca(S_all);

figure;
scatter(score(:,1), score(:,2), 6, '.'); axis equal;





%% Max PC2
[sortedCol, idx] = sort(squeeze(reshape_meaned_score_data(2,:,:)));



FRs = [];

for m = 1:220
    k = idx(m);
    [~,best_fit_index] = min(squeeze(losses(100,2,k,:)));
    FRs = [FRs,sum(sum(squeeze(output(k,best_fit_index,:,:,:))))];
end

%%

%Plotting
for n = 7
    figure('position',[0,0,500,500]);
    
    
    load('all_units_info_with_polished_criteria_modified_perf.mat','all_data');
    load('sound_files.mat','sampleRate','target1','target2');  %Sample Rate 195312 Hz




    
    SpikeTimes = all_data(n).ctrl_tar1_timestamps(:,1);
    
    %The stim lasts from 0s to 2.9801 seconds
    %We are going to run the sim from 0 to 2.9801
    %We need to extract all values between this for the spike distnace loss
    %measures.
    %Switch to zeros for spy plot
    picture = zeros(10,29801);
    
    for m = 1:10
        stim_mask = logical((SpikeTimes{m} > 0) .* (SpikeTimes{m} < 2.9801));
        trial_indicies = round(SpikeTimes{m}(stim_mask)*10000);
        %Switch to one for spy plot
        %1 added for zero indexing. Lets say for instance a spike lands at
        %time = 0. In matlab this would be at the first indicy in the
        %picture which is indicy 1.
        picture(m,trial_indicies+1) = 1;
    end
    
    [~,best_fit_index] = min(losses(100,2,n,:));
    cell_best_fit = squeeze(output(n,best_fit_index,:,:,:));

    figure;
    subplot(4,1,1)
    spy(picture);
    subplot(4,1,2)
    spy(cell_best_fit)
    subplot(4,1,3)
    plot(target1)
    subplot(4,1,4)
    plot(histcounts(mod(find(picture'),29801),149)); hold on
    plot(histcounts(mod(find(cell_best_fit'),29801),149)); 
    legend({'data','best_fit'})
    
    %save('picture_fit'+  string(n) + 'contra' +'.mat','picture');
end

