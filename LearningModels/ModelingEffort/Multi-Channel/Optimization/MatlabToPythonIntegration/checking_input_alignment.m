figure;
subplot(4,1,1)
plot(fr_target_on{1}(2501:32301))
xlim([0 29801])
subplot(4,1,2)
spy(data)
subplot(4,1,3)
spy(squeeze(output(1,:,1,:)))
subplot(4,1,4)
plot(song1(50001:646011))
xlim([0 596011])

%%

figure;
subplot(2,1,1)
mean_mse = mean(losses(:,2,:),3);
plot(mean_mse); hold on
plot(movmean(mean_mse,50))


%%
%Normalize by row????
param_normed = squeeze(params(500,:,:))'./max(squeeze(params(500,:,:))');

figure;

heatmap(param_normed)