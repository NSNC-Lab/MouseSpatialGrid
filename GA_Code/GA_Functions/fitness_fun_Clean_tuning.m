% Fitness function
function fitness = fitness_fun_Clean_tuning(Target_grid,Target_fr_grid, varied_params, optionsGA, plot_all)
    
    %Start by just looking at the R to C connections

    %R to C connections
    varied_struct.RtoC1 = varied_params(1);
    varied_struct.RtoC2 = varied_params(2);
    varied_struct.RtoC3 = varied_params(3);
    varied_struct.RtoC4 = varied_params(4);
    
    % Assume this function updates approximate_grid
    study_dir = fullfile(pwd,'run','4-channel-PV-inputs');
    if exist(fullfile(study_dir, 'solve',['IC_spks_on' '.mat'])) == 2
        load(fullfile(study_dir, 'solve',['IC_spks_on' '.mat']),'spks','dt')
    end
    
    averager = 1;

    %Preallocate so that the appoximate grid can be populated
    approximate_grid = zeros(5,4);
    fr_grid = zeros(5,4);

    for k = 1:averager
        SpikingNetwork_withOffset;
    end

    approximate_grid = approximate_grid/averager;
    fr_grid = fr_grid/averager;
   
    %1. Clean performance
        %Goal here is to enforce the clean condition to be around the average
        %in the target grid.
        
        %a. Going to try to be more strict
        perf_diffs = abs(Target_grid(1,:) - approximate_grid(1,:));
        %clean_cond = mean(Target_grid(1,:));
        
        %b. This should balance out the other fitness function a little bit
        %better hopefully.
        fitness1 = sum(perf_diffs);

        %Scale of 1 to 80

    
    %2. Clean firing rate
        
        %a. Get the firing rate for the target.
        Target_fr = Target_fr_grid(1,:);

        %b. Get the approximate grid firing rate
        approx_fr = fr_grid(1,:);
        
        %c. Penalize for absolute difference
        fitness2 = sum(abs(Target_fr-approx_fr));

        %0 to 160 ish

    %3. Final Fitness calculation
    
    %Fitness weights
    %Might have to adjust these. If we want them of equal importance
    %though, I think this is a good guess to where they woudl be at around
    %equal value given the fitness space is equally complicated for both.
    fw1 = 1;
    fw2 = 8;


    fitness = fitness1*fw1 + fitness2*fw2;

    %disp(['fitness: 1', num2str(fitness1),'  fitness: 2', num2str(fitness2),'   fitness: 3', num2str(fitness3),'   fitness: 4', num2str(fitness4)])
    
    store_values(fitness, approximate_grid,fr_grid, varied_params,optionsGA);

end
