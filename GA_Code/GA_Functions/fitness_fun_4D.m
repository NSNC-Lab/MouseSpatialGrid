% Fitness function
function fitness = fitness_fun_4D(Target_grid,Target_fr_grid, varied_params, optionsGA, plot_all, toggle_real)
    
    %Start by just looking at the R to C connections

    %R to C connections
    varied_struct.RtoC1 = varied_params(1);
    varied_struct.RtoC2 = varied_params(2);
    varied_struct.RtoC3 = varied_params(3);
    varied_struct.RtoC4 = varied_params(4);





    % %REMOVE BEFORE 5D search
    % varied_struct.RtoC1 = 0.00175685264443262;
    % varied_struct.RtoC2 = 0.000198679983008015;
    % varied_struct.RtoC3 = 0.00229808130178778;
    % varied_struct.RtoC4 = 0.00149230691029864;

    
%0.00175685264443262	0.000198679983008015	0.00229808130178778	0.00149230691029864


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
        %tic;
        %profile on;
        SpikingNetwork_withOffset;
        %profile off;
        %profile viewer;
        %toc;
    end

    approximate_grid = approximate_grid/averager;
    fr_grid = fr_grid/averager;
   
    %1. Shape of the tuning curves

    %Need to write in some divide by zero error here, because if the firing
    %rate ends up being zero you get fr_grid_clean_norm = inf and then the
    %output for fitness1 ends up being NaN

    fr_grid_clean_norm = fr_grid(1,:)/(max(fr_grid(1,:)));
    fr_grid_mixed_norm = fr_grid(2,:)/(max(fr_grid(2,:)));

    t_fr_grid_clean_norm = Target_fr_grid(1,:)/max(Target_fr_grid(1,:));
    t_fr_grid_mixed_norm = Target_fr_grid(2,:)/max(Target_fr_grid(2,:));
    
    fitness1 = sum(abs(fr_grid_clean_norm-t_fr_grid_clean_norm)) + sum(abs(fr_grid_mixed_norm-t_fr_grid_mixed_norm));
    
    %Do not allow 0 firing rate to throw a nan
    if(min(fr_grid(2,:)) < 2 || min(fr_grid(1,:)) < 2)
        fitness1 = 1000;
    end
    
    %2. Clean firing rate
    
    fitness2 = 0;
        
    if fr_grid(1,1) > 30
        fitness2 = fitness2 + (fr_grid(1,1) - 30)*20;
    end

    if fr_grid(2,2) < 7
        %For Hz = 6 and below punish by 20 fitness per hz (in lowest area)
        fitness2 = fitness2 + abs(fr_grid(2,4)-7)*20;
    end
    
   %3. Hotspot absolute difference

       fitness3 = abs(Target_grid(2,1) - approximate_grid(2,1));


    %4. Final Fitness calculation
    
    %Fitness weights
    %Might have to adjust these. If we want them of equal importance
    %though, I think this is a good guess to where they woudl be at around
    %equal value given the fitness space is equally complicated for both.
    fw1 = 1;
    fw2 = 1;
    fw3 = 1;

    fitness = fitness1*fw1 + fitness2*fw2 + fitness3*fw3;
   

    %disp(['fitness: 1', num2str(fitness1),'  fitness: 2', num2str(fitness2),'   fitness: 3', num2str(fitness3),'   fitness: 4', num2str(fitness4)])
    
    store_values(fitness, approximate_grid,fr_grid, varied_params,optionsGA,data);

end
