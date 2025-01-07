% Sample data matrix (rows are observations, columns are variables)
data = input_holder;

% Centering the data (optional, as pca function centers automatically)
% data = data - mean(data);

% Perform PCA
[coeff, score, latent, tsquared, explained, mu] = pca(data);

% coeff: principal component coefficients (eigenvectors)
% score: principal component scores (the data projected onto principal components)
% latent: eigenvalues of the covariance matrix (variance explained by each PC)
% explained: percentage of total variance explained by each principal component
% mu: mean of the data used for centering

% Display results
disp('Principal component coefficients (eigenvectors):');
disp(coeff);

disp('Principal component scores:');
disp(score);

disp('Explained variance by each component:');
disp(explained);


figure;
scatter(score(:,1), score(:,2));
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
title('PCA Result');