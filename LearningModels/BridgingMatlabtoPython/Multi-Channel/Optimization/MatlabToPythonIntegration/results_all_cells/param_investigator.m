colors = {'b','g','r','k'};

for k = 1:4

    params_subset = squeeze(params(:,k,:));
    
    %plot(params_subset,colors{k}); hold on

    figure;
    plot(params_subset,colors{k});

end

param_size = size(params);
num_epochs = param_size(1);

%%

figure;
plot(squeeze(sum(params,2)))


%%
a = reshape(params,num_epochs*100,4);
[coeff,score,latent,tsquared,explained] = pca(a);

b = reshape(score,num_epochs,100,4);
figure;
plot(b(:,:,1),'k'); 
figure;
plot(b(:,:,2),'r'); 

%fp1 -> e->e
%fp2 -> e->i
%fp3 -> e->i
%fp4 -> i->e


%quantification of inhibition  fp1 - (fp2->4)

%%

quant_inhib = squeeze(  params(:,1,:) -   (params(:,2,:)+params(:,3,:)+params(:,4,:))/3  );

figure;
plot(quant_inhib)

[c,d] = find(quant_inhib(num_epochs,:)>0);
[e,f] = find(quant_inhib(num_epochs,:)<0);

mostly_inhib = mean(quant_inhib(:,d),2);
mostly_non_inhib = mean(quant_inhib(:,f),2);

figure;
plot(mostly_inhib); hold on
plot(mostly_non_inhib); hold on


%%
params_inhib = params(num_epochs,:,d);
params_non_inhib = params(num_epochs,:,f);


stacked_params = squeeze(cat(3,params_inhib,params_non_inhib));

figure;
heatmap(stacked_params')

avg_params_inhib = squeeze(mean(params_inhib,3));
avg_params_non_inhib = squeeze(mean(params_non_inhib,3));

avg_loss_inhib = squeeze(mean(losses(num_epochs,2,d)));
avg_loss_non_inhib = squeeze(mean(losses(num_epochs,2,f)));

avg_loss_inhib_L2 = squeeze(mean(losses(num_epochs,1,d)));
avg_loss_non_inhib_L2 = squeeze(mean(losses(num_epochs,1,f)));


%%

%example from inhib pop
for k = 77

    figure;
    spy(squeeze(output(k,:,1,:)))

end
%figure;
%spy(squeeze(output(27,:,1,:)))

%%
addpath('C:\Users\ipboy\Documents\GitHub\ModelingEffort\Multi-Channel\Plotting\OliverDataPlotting')
addpath('C:\Users\ipboy\Documents\GitHub\ModelingEffort\Multi-Channel\Plotting\OliverDataPlotting\PicturesToFit')
load('all_units_info_with_polished_criteria_modified_perf.mat','all_data');
data = load('picture_fit' + string(7) + 'contra.mat').picture;

figure;
subplot(3,1,1)
spy(data)
title('Data (Ground truth)')
subplot(3,1,2)
spy(squeeze(output(77,:,1,:)))
title('FR: 2.5Hz -- MSE PSTH Loss :  ' +string(losses(num_epochs,2,77)))
subplot(3,1,3)
spy(squeeze(output(32,:,1,:)))
title('FR: 32Hz -- MSE PSTH Loss :  ' +string(losses(num_epochs,2,32)))

FRS = sum(sum(squeeze(output),2),3);

[idx,val] = min(abs(608-FRS))

%%
figure;
%plot(squeeze(params))


figure;
mse = squeeze(losses(:,2,:))';

plot(mean(mse)); hold on
plot(movmean(mean(mse),100)); hold on

figure;
plot(movmean(mean(mse),100));

%% Pearsons correlation
num_epochs = 100;
loss_vals = (squeeze(losses(num_epochs,2,:)));
last_params = (squeeze(params(num_epochs,:,:)))';

%Trying across all trials
%reshape_loss_vals = (reshape(losses(:,2,:),[num_epochs*100,1]));
%reshape_params = (reshape(params,[num_epochs*100,13]));

%R = corr([reshape_params reshape_loss_vals]);  

R = corr([last_params loss_vals]);       % correlation matrix
corrs = R(1:end-1,end)

figure;
heatmap(corrs)
n = 256;                % number of colors
half = n/2;

g = [linspace(0,1,half), linspace(1,0,half)]';  % dark → light → dark
cmap = [g g g];                                  % make it grayscale

colormap(cmap);
colorbar;

%%
%names = {'gsyn On->ROn','gsyn On->SOnOff','gsyn Off->SOnOff','gsyn SOnOff->ROn','FR','tau ad On','tau ad Off','tau ad SOnOff','tau ad ROn','g inc On','g inc Off','g inc SOnOff','g inc ROn'}; 
names = {'gsyn On->ROn','gsyn On->SOnOff','gsyn Off->SOnOff','gsyn SOnOff->ROn','FR','tau ad On','tau ad Off','tau ad SOnOff','tau ad ROn','g inc On','g inc Off','g inc SOnOff','g inc ROn','tauP On','tauP Off','tauP SOnOff','TauP Ron'}; 

figure
for k = 1:17
    subplot(4,5,k)
    plot(squeeze(params(:,k,:)),'k')
    title(names{k})
end


%% PCA of output

a = squeeze(params(num_epochs,:,:))'./max(squeeze(params(num_epochs,:,:))');

[coeff,score,latent,tsquared,explained] = pca(a);


figure;
heatmap(coeff)

figure;
scatter(score(:,1),score(:,2), 100, squeeze(losses(num_epochs,2,:)), 'filled')
colormap(parula);   % or your custom colormap
colorbar;

xlabel('PC1 (As this gets more positive - Onset Adaption increases and Relay adaption decreases)')
ylabel('PC2 (As this gets more positive - Offset Adaption increases and PV adaption decreases)')


%% Parameter movement
paramsa = reshape(params,[num_epochs 100 13]);
% param_hist: [E x T x P]
[E, T, P] = size(paramsa);

final_params = squeeze(paramsa(end,:,:));  % [T x P]
d = zeros(E,1);

for e = 1:E
    params_e = squeeze(paramsa(e,:,:));   % [T x P]
    diff = params_e - final_params;          % [T x P]
    d(e) = mean( sqrt(sum(diff.^2, 2)) );    % mean L2 over trials
end

d_norm = d / d(1);

figure;
plot(1:E, d_norm, '-o');
xlabel('Epoch');
ylabel('Normalized distance to final params');
title('Convergence of parameters (distance to final epoch)');
grid on;

%Per parameter
[E, T, P] = size(paramsa);

final_params = squeeze(paramsa(end, :, :));   % [T x P]
dist_to_final = zeros(E, P);                     % [E x P]

for e = 1:E
    params_e = squeeze(paramsa(e, :, :));     % [T x P]
    diff = params_e - final_params;              % [T x P]
    
    % mean |difference| across trials for each parameter
    dist_to_final(e, :) = mean(abs(diff), 1);    % 1 x P
end

% Normalize per parameter so epoch 1 = 1
dist_norm = dist_to_final ./ (dist_to_final(1, :) + eps);

% Subplot layout
nrows = ceil(sqrt(P));
ncols = ceil(P / nrows);

figure;
for p = 1:P
    subplot(nrows, ncols, p);
    plot(1:E, dist_norm(:, p), 'LineWidth', 1.5);
    xlabel('Epoch');
    ylabel('Norm dist');
    title(sprintf('Param %d', p));
    grid on;
end
sgtitle('Distance to final value per parameter');


%% Parameter Step size
step = zeros(E-1,1);
eps_val = 1e-8;

for e = 2:E
    params_prev = squeeze(paramsa(e-1,:,:));   % [T x P]
    params_cur  = squeeze(paramsa(e,:,:));     % [T x P]

    diff = params_cur - params_prev;              % [T x P]
    num  = sqrt(sum(diff.^2, 2));                 % [T x 1]
    denom = sqrt(sum(params_prev.^2, 2)) + eps_val;

    step(e-1) = mean(num ./ denom);               % average over trials
end

figure;
plot(2:E, step, '-o');
xlabel('Epoch');
ylabel('Relative step size');
title('Parameter update size per epoch');
grid on;

step_param = zeros(E-1, P);
eps_val = 1e-8;

for e = 2:E
    params_prev = squeeze(paramsa(e-1, :, :));  % [T x P]
    params_cur  = squeeze(paramsa(e,   :, :));  % [T x P]
    
    diff  = params_cur - params_prev;              % [T x P]
    
    % mean |delta| across trials
    num   = mean(abs(diff), 1);                    % 1 x P
    % mean |prev| across trials (for relative size)
    denom = mean(abs(params_prev), 1) + eps_val;   % 1 x P
    
    step_param(e-1, :) = num ./ denom;             % relative step per param
end

% Same subplot layout as above
figure;
for p = 1:P
    subplot(nrows, ncols, p);
    plot(2:E, step_param(:, p), 'LineWidth', 1.5);
    xlabel('Epoch');
    ylabel('Rel step');
    title(sprintf('Param %d', p));
    grid on;
end
sgtitle('Relative step size per parameter');

%%
 figure;
subplot(2,1,1)
spy(data)
subplot(2,1,2)
spy(squeeze(best_output))


%%
figure;
plot(squeeze(min(losses(:,2,:),[],3))); hold on
plot(movmean(squeeze(min(losses(:,2,:),[],3)),100))

figure;
plot(movmean(squeeze(min(losses(:,2,:),[],3)),100))