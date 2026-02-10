figure;

colors = {'b','g','r','k','c'};

%params_no_NAN = params(:,:,1:97);

for k = 1:5
    plot(mean(params(:,k,:),3),colors{k}); hold on
    %plot(mean(params_no_NAN(:,k,:),3),colors{k}); hold on
end
legend({'On->ROn','Off->ROn','On->PV','Off->PV','PV->ROn'})
%%
figure;
spy(best_output)

%%
figure;
spy(squeeze(output(1,:,:,:)))

%%
figure
subplot(3,1,1)
plot(min(squeeze(losses(4:end,2,:))'))
title('min-PSTH-loss over epochs')
subplot(3,1,2)
plot(mean(squeeze(losses(:,2,:))'))
title('avg-PSTH-loss over epochs')
subplot(3,1,3)
plot(squeeze(mean(sum(abs(diff(params)),2),3)))
%plot(squeeze(mean(mean(diff(params_no_NAN),2),3)))
title('convergence')

%%
cd(userpath);
cd('../GitHub/ModelingEffort/Multi-Channel/Plotting/OliverDataPlotting')

load('all_units_info_with_polished_criteria_modified_perf.mat','all_data');
load('sound_files.mat','sampleRate','target1','target2');  %Sample Rate 195312 Hz


cd('PicturesToFit\')
%Look at the first animal type
n = 7;
    
SpikeTimes = all_data(n).ctrl_tar1_timestamps(:,1);

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


figure;
spy(picture)


%%
figure
spy(data)

%%
figure
spy(forwards_out)

%%
figure
spy(squeeze(noise(1,:,:)))

%%
figure;
plot(histcounts(find(output(:,:,:)),149),'b','LineWidth',2); hold on

%%
plot(histcounts(find(output(:,:,:)),149),'g','LineWidth',2); hold on

%%
plot(histcounts(find(data(:,:)),149)*100,'k','LineWidth',2); hold on
%%
legend('nsyn=0.011','nsyn=0.015','data')


%%
figure;
spy(squeeze(noise(3,:,:)))

%%
mean(losses(100,2,:))

%%
X = squeeze(params(100,:,:))';
Xz = zscore(X);
[coeff, score] = pca(Xz);

pc = score(:,1:2);

% Grid for density evaluation
xg = linspace(min(pc(:,1)), max(pc(:,1)), 200);
yg = linspace(min(pc(:,2)), max(pc(:,2)), 200);
[Xg, Yg] = meshgrid(xg, yg);

% 2D kernel density estimate on grid
F = ksdensity(pc, [Xg(:), Yg(:)]);
F = reshape(F, size(Xg));

figure; hold on;
scatter(pc(:,1), pc(:,2), 200, 'k.');   % your scatter
%contour(Xg, Yg, F, 10, 'LineWidth', 1);   % topology lines
contourf(Xg, Yg, F, 10, 'LineStyle', 'none','FaceAlpha',0.3);
% or: contourf(Xg, Yg, F, 10, 'LineStyle', 'none'); alpha(0.3);
hold off;
xlabel('PC1'); ylabel('PC2');

%%
E = size(params,1); P = size(params,2); B = size(params,3);
t = 1:E;
Xfinal = squeeze(params(end,:,:))';     % [B x P]
X = params;                             % shorthand

mu = squeeze(mean(X,3));              % [E x P]
sd = squeeze(std(X,0,3));             % [E x P]

figure;
for p=1:P
    subplot(P,1,p); hold on;
    fill([t fliplr(t)], [mu(:,p)'+sd(:,p)' fliplr(mu(:,p)'-sd(:,p)')], ...
         'k', 'FaceAlpha',0.1, 'EdgeColor','none');
    plot(t, mu(:,p), 'k', 'LineWidth',1.5);
    ylabel(sprintf('p%d',p));
    if p==1, title('Mean \pm SD across batch'); end
    if p==P, xlabel('epoch'); end
    hold off;

end

%%
E = size(params,1); P = size(params,2); B = size(params,3);
t = 1:E;
Xfinal = squeeze(params(end,:,:))';     % [B x P]
X = params;                             % shorthand

mu = squeeze(mean(X,3));              % [E x P]
sd = squeeze(std(X,0,3));             % [E x P]

figure;
for p=1:P
    %subplot(P,1,p); hold on;
    fill([t fliplr(t)], [mu(:,p)'+sd(:,p)' fliplr(mu(:,p)'-sd(:,p)')], ...
         colors{p}, 'FaceAlpha',0.1, 'EdgeColor','none'); hold on;
    plot(t, mu(:,p), colors{p}, 'LineWidth',1.5); hold on;
    ylabel(sprintf('p%d',p));
    if p==1, title('Mean \pm SD across batch'); end
    if p==P, xlabel('epoch'); end
    %hold off;

end

%%

q10 = squeeze(prctile(X,10,3));
q50 = squeeze(prctile(X,50,3));
q90 = squeeze(prctile(X,90,3));

figure;
for p=1:P
    subplot(P,1,p); hold on;
    fill([t fliplr(t)], [q90(:,p)' fliplr(q10(:,p)')], 'k', ...
         'FaceAlpha',0.1, 'EdgeColor','none');
    plot(t, q50(:,p), 'k', 'LineWidth',1.5);
    ylabel(sprintf('p%d',p));
    if p==1, title('Median with 10–90% band'); end
    if p==P, xlabel('epoch'); end
    hold off;
end

%%
nShow = 20;
idx = randperm(B, nShow);

figure;
for p=1:P
    subplot(P,1,p); hold on;
    for k=1:nShow
        plot(t, squeeze(X(:,p,idx(k))), 'Color',[0 0 0 0.15]); % light lines
    end
    plot(t, squeeze(mean(X(:,p,:),3)), 'k', 'LineWidth',1.5);
    ylabel(sprintf('p%d',p));
    if p==1, title('Batch trajectories (subset) + mean'); end
    if p==P, xlabel('epoch'); end
    hold off;
end

%%
p = 3; % choose param
M = squeeze(X(:,p,:))';   % [B x E]

figure;
imagesc(t, 1:B, M); axis tight;
xlabel('epoch'); ylabel('batch');
title(sprintf('Param %d over training (batch x epoch)', p));
colorbar;


%%
C = corr(Xfinal);   % [P x P]
figure;
imagesc(C); axis square;
xticks(1:P); yticks(1:P);
title('Correlation of parameters (final epoch)');
colorbar;

%%
figure;
plotmatrix(Xfinal);
title('Final epoch: pairwise scatter of parameters');

%%
p1=1; p2=2;
pc2 = [Xfinal(:,p1), Xfinal(:,p2)];

xg = linspace(min(pc2(:,1)), max(pc2(:,1)), 200);
yg = linspace(min(pc2(:,2)), max(pc2(:,2)), 200);
[Xg,Yg] = meshgrid(xg,yg);
F = ksdensity(pc2, [Xg(:),Yg(:)]);
F = reshape(F, size(Xg));

figure; hold on;
scatter(pc2(:,1), pc2(:,2), 25, 'o');
contour(Xg,Yg,F,10,'LineWidth',1);
xlabel(sprintf('p%d',p1)); ylabel(sprintf('p%d',p2));
title('Final epoch: scatter + density contours');
hold off;

%%
Xz = zscore(Xfinal);
[coeff, score, ~, ~, explained] = pca(Xz);

figure;
scatter(score(:,1), score(:,2), 25, 'o');
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title('Final epoch PCA (batch samples)');

%% -- This one looks cool
Xfinal_z = zscore(Xfinal);
[coeff,~,~,~,explained,mu] = pca(Xfinal_z);

S = zeros(B,2,E);
for e=1:E
    Xe = squeeze(X(e,:,:))';      % [B x P]
    Xez = (Xe - mean(Xfinal,1)) ./ std(Xfinal,0,1); % normalize using final stats
    Se = Xez * coeff(:,1:2);
    S(:,:,e) = Se;
end

figure; hold on;
for e=[1 50 200 500 1000]
    scatter(S(:,1,e), S(:,2,e), 15, 'o');
end
xlabel('PC1'); ylabel('PC2');
title('Projected PCA space at selected epochs');
legend('e1','e50','e200','e500','e1000');
hold off;

%%
D = zeros(E,B);
for e=1:E
    Xe = squeeze(X(e,:,:))';          % [B x P]
    D(e,:) = vecnorm(Xe - Xfinal, 2, 2);  % L2 per batch sample
end

figure;
plot(t, mean(D,2), 'k', 'LineWidth',1.5);
xlabel('epoch'); ylabel('mean ||p(e)-p(final)||');
title('Convergence to final parameters (mean across batch)');

%%
N = zeros(E,B);
for e=1:E
    Xe = squeeze(X(e,:,:))';
    N(e,:) = vecnorm(Xe,2,2);
end
figure;
plot(t, mean(N,2), 'k', 'LineWidth',1.5);
xlabel('epoch'); ylabel('mean ||p||');
title('Parameter vector norm over training');


%%
% params: [epochs x params x batch] = [1000 x 5 x 100]

E = size(params,1);
P = size(params,2);

figure;
tl = tiledlayout(1,P,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Final epoch (epoch %d) parameter distributions', E));

for p = 1:P
    x = squeeze(params(end,p,:));   % [batch x 1]

    nexttile; hold on;

    % Histogram as PDF
    histogram(x, 'Normalization','pdf','NumBins',30);

    % Smooth density overlay (KDE)
    [f,xi] = ksdensity(x);
    plot(xi, f, 'LineWidth', 1.5);

    % Median marker
    med = median(x);
    yl = ylim;
    plot([med med], yl, 'k--', 'LineWidth', 1);

    title(sprintf('Param %d', p));
    xlabel('value'); ylabel('pdf');
    hold off;
end
%%

% params: [epochs x params x batch] = [E x P x B]
E = size(params,1);
P = size(params,2);
B = size(params,3);

Xfinal = squeeze(params(end,:,:))';   % [B x P]  (rows=batch, cols=params)

% Sort batch samples by parameter #3 (ascending)
[~, idx] = sort(Xfinal(:,3), 'ascend');
Xs = Xfinal(idx,:);                  % [B x P] sorted

% --- Heatmap (raw values) ---
figure;
imagesc(Xs');                        % [P x B]
axis tight;
xlabel('batch (sorted by param 3)');
ylabel('parameter');
yticks(1:P);
yticklabels("p" + (1:P));
title(sprintf('Final epoch heatmap (epoch %d), sorted by param 3', E));
colorbar;

% --- Optional: heatmap with per-parameter z-scoring (often more readable) ---
% Xz = (Xs - mean(Xs,1)) ./ std(Xs,0,1);  % normalize each parameter across batch
% figure;
% imagesc(Xz');
% axis tight;
% xlabel('batch (sorted by param 3)');
% ylabel('parameter (z-scored)');
% yticks(1:P); yticklabels("p" + (1:P));
% title(sprintf('Final epoch heatmap (z-scored), sorted by param 3', E));
% colorbar;

%%
% params: [E x P x B]
% loss:   [E2 x B]   (use loss(end,:) as the final-epoch loss per batch member)

E = size(params,1);
P = size(params,2);
B = size(params,3);

loss_final = squeeze(losses(end,2,:));            % [B x 1]

nbins = 10;                            % granularity (change this)

figure;
tl = tiledlayout(1,P,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Final epoch (epoch %d): hist + mean loss per bin', E));

for p = 1:P
    x = squeeze(params(end,p,:));      % [B x 1]

    % Define bins (you can also use 'BinWidth' or 'BinMethod')
    edges = linspace(min(x), max(x), nbins+1);
    bin   = discretize(x, edges);      % [B x 1], bin index in 1:nbins

    % Mean loss per bin
    fillv = cast(NaN, 'like', loss_final);
    meanLoss = accumarray(bin(valid), loss_final(valid), [nbins 1], @mean, fillv);

    % (Optional) SEM per bin
    %nPerBin = accumarray(bin(valid), 1, [nbins 1], @sum, 0);
    %stdLoss = accumarray(bin(valid), loss_final(valid), [nbins 1], @std, NaN);
    %semLoss = stdLoss ./ sqrt(max(nPerBin,1));

    centers = (edges(1:end-1) + edges(2:end)) / 2;

    nexttile; hold on;

    % Histogram (left y-axis)
    yyaxis left;
    histogram(x, 'BinEdges', edges, 'Normalization','pdf');
    ylabel('pdf');

    % Mean loss curve (right y-axis)
    yyaxis right;
    plot(centers, meanLoss, 'k-', 'LineWidth', 1.5);
    % Optional errorbars:
    % errorbar(centers, meanLoss, semLoss, 'k.', 'CapSize', 0);

    ylabel('mean loss');

    xlabel('value');
    title(sprintf('Param %d', p));
    hold off;
end