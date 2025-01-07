function [counts, spont, PSTH1] = calc_num_spikes(timeStamps1,nbins)

 concatenatedTimestamps = cell(1, 4);
    counts = zeros(size(timeStamps1));
    spont = zeros(size(timeStamps1));
    PSTH1 = [];
    % Loop through each cell to count elements in the range [0, 3] and [-1 0]
    for j = 1:size(timeStamps1, 2) % Loop over columns
    
        tempTimestamps = [];
        for i = 1:size(timeStamps1, 1) % Loop over rows
            % Extract the data from the current cell
            currentData = timeStamps1{i, j};
             tempTimestamps = [tempTimestamps;currentData];
            % Count the number of elements in the range [0, 3]
            counts(i, j) = sum(currentData >= 0 & currentData <= 3);
            spont(i,j) = sum(currentData >= -1 & currentData < -0.5) + sum (currentData >= 3.5);
        end
        concatenatedTimestamps{j} = sort(tempTimestamps);

        [PSTH1(j,:),edges] = histcounts(concatenatedTimestamps{j},nbins); 
    end

end