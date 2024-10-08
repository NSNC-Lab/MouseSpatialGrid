grid_tracker = {};
fr_tracker = {};

epochs = 200;

for k = 1:epochs
    SpikingNetwork_withOffset;
    grid_tracker{end+1} = approximate_grid;
    fr_tracker{end+1} = fr_grid;
    k/epochs
end

%Look at the scenario where inhibition is not as intense as usual