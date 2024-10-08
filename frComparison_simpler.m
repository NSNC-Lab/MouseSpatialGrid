%Find Associated Ron FRs in the output
arrays = cellfun(@(x) x.ROn, {data.fr}, 'UniformOutput', false);

% Concatenate all the arrays into a 4x24 array
flattenedArray = vertcat(arrays{:});
output = [flattenedArray.channel1;flattenedArray.channel2;flattenedArray.channel3;flattenedArray.channel4];

save("Baseline.mat","output")