%%
figure;
plot(squeeze(losses(:,2,:)))

%% Perform PCA

%params_transposed = reshape(params,[64,16384,5]);
%params_reshaped = reshape(permute(params,[2 1 3]), [], 5);

params_perm1 = permute(params,[3,1,2]);
params_reshaped = reshape(params_perm1,1048576,5);
params_reshaped_2 = reshape(params_reshaped,[16,16,16,16,16,5]);

losses_perm1 = permute(losses,[3,1,2]);
losses_reshaped = reshape(losses_perm1,1048576,2);
losses_reshaped_2 = reshape(losses_reshaped,[16,16,16,16,16,2]);
lossses_perm2 = permute(losses_reshaped_2,[6,1,2,3,4,5]);
lossses_mse = lossses_perm2(2,:,:,:,:,:);

J1 = squeeze(min(lossses_mse, [], [2 3 4]));   % size 16x1
%plot(p1_vals, J1, '-o'); xlabel('p1'); ylabel('best score')



%Xfinal_z = zscore(params_reshaped);
%[coeff,~,~,~,explained,mu] = pca(Xfinal_z);


%%
figure;
plot(losses_reshaped(:,1))
figure;
plot(losses_reshaped(:,2))


%%
