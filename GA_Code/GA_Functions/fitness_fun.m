% Fitness function
function fitness = fitness_fun(Target_grid,Target_fr_grid, varied_params, optionsGA)
    
    %Intra Channel PV connections
    varied_struct.IntraPV1 = varied_params(1);
    %varied_struct.IntraPV2 = varied_params(2);
    %varied_struct.IntraPV3 = varied_params(3);
    %varied_struct.IntraPV4 = varied_params(4);

    %Outer Cross cannel PV connections
    varied_struct.CrossPV1 = varied_params(2);
    varied_struct.CrossPV2 = varied_params(3);
    varied_struct.CrossPV3 = varied_params(4);
    %varied_struct.CrossPV4 = varied_params(8);
    %varied_struct.CrossPV5 = varied_params(9);
    %varied_struct.CrossPV6 = varied_params(10);
    %varied_struct.CrossPV7 = varied_params(11);
    %varied_struct.CrossPV8 = varied_params(12);
    %varied_struct.CrossPV9 = varied_params(13);
    %varied_struct.CrossPV10 = varied_params(14);
    %varied_struct.CrossPV11 = varied_params(15);
    %varied_struct.CrossPV12 = varied_params(16);

    %On to R connections
    %varied_struct.ONtoR1 = varied_params(17);
    %varied_struct.ONtoR2 = varied_params(18);
    %varied_struct.ONtoR3 = varied_params(19);
    %varied_struct.ONtoR4 = varied_params(20);

    %R to C connections
    varied_struct.RtoC1 = varied_params(5);
    %varied_struct.RtoC2 = varied_params(22);
    %varied_struct.RtoC3 = varied_params(23);
    varied_struct.RtoC4 = varied_params(6);

    %Cross channel SOM connections (From Contra)
    %varied_struct.XR1 = varied_params(25);
    %varied_struct.XR2 = varied_params(26);
    %varied_struct.XR3 = varied_params(27);
    %varied_struct.XR4 = varied_params(28);
    %varied_struct.XR5 = varied_params(29);
    %varied_struct.XR6 = varied_params(30);
    %varied_struct.XR7 = varied_params(31);
    %varied_struct.XR8 = varied_params(32);
    %varied_struct.XR9 = varied_params(33);
    %varied_struct.XR10 = varied_params(34);
    %varied_struct.XR11 = varied_params(35);
    %varied_struct.XR12 = varied_params(36);
    
    % Assume this function updates approximate_grid
    study_dir = fullfile(pwd,'run','4-channel-PV-inputs');
    if exist(fullfile(study_dir, 'solve',['IC_spks_on' '.mat'])) == 2
        load(fullfile(study_dir, 'solve',['IC_spks_on' '.mat']),'spks','dt')
    end
    

    averager = 1;
    %Fr_Frac = 0.3;
    %Perf_Frac = 1;
    
    %Preallocate so that the appoximate grid can be populated
    approximate_grid = zeros(5,4);
    fr_grid = zeros(5,4);

    for k = 1:averager
        SpikingNetwork_withOffset;
    end

    approximate_grid = approximate_grid/averager;
    fr_grid = fr_grid/averager;
    
    %Putting firing rate and performance together.
    %fitness = sum(sum(abs(Target_fr_grid(1:2,:)-fr_grid(1:2,:))))*Fr_Frac...
    %    + abs(Target_grid(2,1) - (approximate_grid(2,1)))*Perf_Frac;
    %fitness = abs(Target_grid(2,1) - (approximate_grid(2,1))); Perf
    %fitness = sum(sum(abs(Target_fr_grid(1:2,:)-fr_grid(1:2,:)))) FR
    

    %New More relaxed Fitness 8/6
    
    %1. Mixed condition
    
        %Perfect score is when the upper left quadrent is clearly above
        %baseline
        
        %a. Find background
    
        background = mean(mean(Target_grid(2:5,:)));
    
        %find the points above or below backgorund the upper left quadrent is
        %This will have to be changed if we are not looking at the upper left
        %quadrent.
        %Perhaps there is some way to loclize this so taht it always does an
        %average around the hotspot.
        
        %b find hotspot
        %8/8 just changed to be the upper left for now so that it is more
        %focused
        %hotspot = mean(mean(Target_grid(2:3,1:2)));
        hotspot = mean(mean(Target_grid(2,1)));
    
        %c find level above background we are trying to recreate
    
        fit_level = hotspot-background;
        
        %d. Subtract our estimate from the fit_level
        
        background_approx = mean(mean(approximate_grid(2:5,:)));
        hotspot_approx = mean(mean(approximate_grid(2:3,1:2)));
    
        fit_level_approximate = hotspot_approx-background_approx;
    
        %e. Penelize it the further it is away.

        fitness1 = (fit_level-fit_level_approximate)*0.1;

    %2. Firing Rate Match

        %Here we are just making sure that the firing rate is constrained
        %Goal here is to just be less than 30 Hz
        %Penalize things over this value.
        fitness2 = sum(fr_grid(fr_grid>30) - 30)*0.1;

    %3. Relative Firing rate match

        %a. Take the clean condition firing rates
        %Compute the ratio of the difference between the two curves
        Target_fr_diff = Target_fr_grid(1,:)./Target_fr_grid(2,:);

        %b. Penalize the differential between the appoximate and the target
        Resultant_diff = fr_grid(1,:)./fr_grid(2,:);
        
        %Used to be sum() but you run into an isse when mixed condition
        %goes to 0. Still think this is important though, it might need
        %large weight in future. Added abs value so that getting to far
        %below the ratio required causes a worse fit.
        fitness3 = sum(abs(Target_fr_diff-Resultant_diff));

    %4. Clean condition for grid

        %a. Get the clean firing rates
        clean_cond = mean(Target_grid(1,:));

        %b. Penalize them from being away from the target
        fitness4 = (clean_cond - mean(approximate_grid(1,:)))*0.05;

    %5. Final Fitness calculation
    
    %Fitness weights
    fw1 = 5;
    fw2 = 0.5;
    fw3 = 2;
    fw4 = 1;


    fitness = fitness1*fw1 + fitness2*fw2 + fitness3*fw3 + fitness4*fw4;

    %disp(['fitness: 1', num2str(fitness1),'  fitness: 2', num2str(fitness2),'   fitness: 3', num2str(fitness3),'   fitness: 4', num2str(fitness4)])
    
    store_values(fitness, approximate_grid,fr_grid, varied_params,optionsGA);

end
