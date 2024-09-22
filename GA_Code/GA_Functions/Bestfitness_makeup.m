% bestApproximateGrid = bestApproximateGrid2;
% bestfrGrid = bestfrGrid2;


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
    hotspot = mean(mean(Target_grid(2:3,1:2)));
    

    %c find level above background we are trying to recreate

    fit_level = hotspot-background;
    
    %d. Subtract our estimate from the fit_level
    
    background_approx = mean(mean(bestApproximateGrid(2:5,:)));
    hotspot_approx = mean(mean(bestApproximateGrid(2:3,1:2)));

    fit_level_approximate = hotspot_approx-background_approx;

    %e. Penelize it the further it is away.

    fitness1 = (fit_level-fit_level_approximate)*0.1;

%2. Firing Rate Match

    %Here we are just making sure that the firing rate is constrained
    %Goal here is to just be less than 30 Hz
    %Penalize things over this value.
    fitness2 = sum(bestfrGrid(bestfrGrid>30) - 30)*0.1;

%3. Relative Firing rate match

    %a. Take the clean condition firing rates
    %Compute the ratio of the difference between the two curves
    Target_fr_diff = Target_fr_grid(1,:)./Target_fr_grid(2,:);

    %b. Penalize the differential between the appoximate and the target
    Resultant_diff = bestfrGrid(1,:)./bestfrGrid(2,:);
    
    %Used to be sum() but you run into an isse when mixed condition
    %goes to 0. Still think this is important though, it might need
    %large weight in future. Added abs value so that getting to far
    %below the ratio required causes a worse fit.
    fitness3 = sum(abs(Target_fr_diff-Resultant_diff));

%4. Clean condition for grid

    %a. Get the clean firing rates
    clean_cond = mean(Target_grid(1,:));

    %b. Penalize them from being away from the target
    fitness4 = (clean_cond - mean(bestApproximateGrid(1,:)))*0.05;

 %5. Final Fitness calculation

fitness = fitness1 + fitness2 + fitness3 + fitness4;
