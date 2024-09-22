% bestApproximateGrid = bestApproximateGrid2;
% bestfrGrid = bestfrGrid2;


%1. Clean performance
    %Goal here is to enforce the clean condition to be around the average
    %in the target grid.
    
    %a. Get the mean firing performance values for the target
    clean_cond = mean(Target_grid(1,:));
    
    %b. Penalize them from being away from the target
    %Going to put an abs in for now because we know that for the easy
    %condtion we can get way over the values that we are after. This
    % will make it so that we try to get our grids a little bit better.
    fitness1 = abs(clean_cond - mean(bestApproximateGrid(1,:)));

    %Scale of 0 to 20


%2. Clean firing rate
    
    %a. Get the firing rate for the target.
    Target_fr = Target_fr_grid(1,:);

    %b. Get the approximate grid firing rate
    approx_fr = bestfrGrid(1,:);
    
    %c. Penalize for absolute difference
    fitness2 = sum(abs(Target_fr-approx_fr));

    %0 to 160 ish


%3. Hotspot absolute difference

    fitness3 = abs(Target_grid(2,1) - bestApproximateGrid(2,1));


%4. Final Fitness calculation

%Fitness weights
%Might have to adjust these. If we want them of equal importance
%though, I think this is a good guess to where they woudl be at around
%equal value given the fitness space is equally complicated for both.
fw1 = 1;
fw2 = 8;
fw3 = 8;


fitness = fitness1*fw1 + fitness2*fw2 + fitness3*fw3;
