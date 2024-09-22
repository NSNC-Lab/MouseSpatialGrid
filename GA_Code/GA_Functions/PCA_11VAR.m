%Perform PCA
% Standardize the data
[data_standardized, mu, sigma] = zscore(Total_search_fitness);

% Perform PCA
[coeff, score, latent] = pca(data_standardized);



%Here we can do a couple of different things. Perhaps lets try plotting the
%evlolution of the mean variable place on the first 2 PCs.

PCA_evolution = [];

for c = 1:nVars
    cur_gen = Total_search_fitness(c,:);
    cur_gen_standardized = (cur_gen - mu) ./ sigma;
    cur_gen_projection_mean = max(cur_gen_standardized * transpose(coeff(1:2, :)));
    PCA_evolution = [PCA_evolution; cur_gen_projection_mean];
end

figure;
scatter(PCA_evolution(:,1),PCA_evolution(:,2))
