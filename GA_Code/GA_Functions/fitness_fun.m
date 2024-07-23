% Fitness function
function fitness = fitness_fun(Target_grid, varied_params, optionsGA)
    
    varied_struct.IntraPV = varied_params(1);
    varied_struct.RtoC = varied_params(2);
    varied_struct.XR = varied_params(3);
    
    % Assume this function updates approximate_grid
    
    averager = 1;
    
    %Preallocate so that the appoximate grid can be populated
    approximate_grid = zeros(5,4);

    for k = 1:averager
        SpikingNetwork_withOffset;
    end

    approximate_grid = approximate_grid/averager;
   
    fitness = abs(Target_grid(2,1) - (approximate_grid(2,1)));
    
    store_values(fitness, approximate_grid, varied_params,optionsGA);

end
