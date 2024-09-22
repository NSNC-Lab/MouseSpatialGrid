% Fitness function
function fitness = fitness_2D(Target_grid,Target_fr_grid, varied_params, optionsGA, plot_all)
    
    %Start by just looking at the R to C connections

    %R to C connections
    varied_struct.RtoC1 = varied_params(1);
    varied_struct.RtoC4 = varied_params(2);

    
    % Assume this function updates approximate_grid
    study_dir = fullfile(pwd,'run','4-channel-PV-inputs');
    if exist(fullfile(study_dir, 'solve',['IC_spks_on' '.mat'])) == 2
        load(fullfile(study_dir, 'solve',['IC_spks_on' '.mat']),'spks','dt')
    end
    

    %Preallocate so that the appoximate grid can be populated
    approximate_grid = zeros(5,4);
    fr_grid = zeros(5,4);

    SpikingNetwork_withOffset;

   
    %Goal here is to match firing rates
    fitness = abs(Target_fr_grid(1,1)-fr_grid(1,1));

    %disp(['fitness: 1', num2str(fitness1),'  fitness: 2', num2str(fitness2),'   fitness: 3', num2str(fitness3),'   fitness: 4', num2str(fitness4)])
    
    store_values(fitness, approximate_grid,fr_grid, varied_params,optionsGA);

end
