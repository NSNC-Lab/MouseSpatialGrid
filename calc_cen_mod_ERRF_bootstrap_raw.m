function [errf_width_ctrl_all, mod_depth_ctrl_all, centroid_ctrl_all, bootstrapMeans]= calc_cen_mod_ERRF_bootstrap_raw(numBootstraps, counts,degree_list_x)  


    % Initialize the matrix to hold the bootstrapped column-wise averages
    bootstrapMeans = zeros(numBootstraps, 4);
    
    % Perform bootstrapping
    for i = 1:numBootstraps
        % Sample with replacement from each column
        sampledMatrix = datasample(counts, size(counts, 1), 1, 'Replace', true);
        
        % Calculate the mean of each column in the sampled matrix
        bootstrapMeans(i, :) = mean(sampledMatrix, 1);
        mod_depth_ctrl_all(i,:) = (max(bootstrapMeans(i, :)) - min(bootstrapMeans(i, :)))/max(bootstrapMeans(i, :))*100;
        centroid_ctrl_all(i,:) = calc_centroid_Oliver(degree_list_x, bootstrapMeans(i, :));
        [~,errf_width_ctrl_all(i,:)]= calc_errf(flip(degree_list_x),flip(bootstrapMeans(i, :))); % flip the order, otherwise it will generate a negative number
    end

    