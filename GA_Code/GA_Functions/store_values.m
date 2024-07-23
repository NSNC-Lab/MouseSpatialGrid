% Keep persistent memory of variables in the GA
function store_values(fitness, approximate_grid, varied_params, optionsGA,state)

    persistent bestFitness bestApproximateGrid bestvars curfitness curvars curgrid

    if isempty(bestFitness)
        bestFitness = inf;
        bestApproximateGrid = [];
        bestvars = [];
        curvars = [];
        curfitness = [];
        curgrid = {};
    end

    if(length(curfitness) == optionsGA.PopulationSize)
        curvars = [];
        curfitness = [];
        curgrid = {};
    end

    % Update the best fitness and best approximate grid if current fitness is better
    if fitness < bestFitness
        bestFitness = fitness;
        bestApproximateGrid = approximate_grid;
        bestvars = varied_params;
    end
    

    curvars = [curvars varied_params];
    curfitness = [curfitness fitness];
    curgrid = {curgrid approximate_grid};


    % Save the variables back to the main workspace
    assignin('base', 'bestFitness', bestFitness);
    assignin('base', 'bestApproximateGrid', bestApproximateGrid);
    assignin('base', 'bestvars', bestvars);
    assignin('base', 'curvars', curvars);
    assignin('base', 'curfitness', curfitness);
    assignin('base', 'curgrid', curgrid);

end