% Keep persistent memory of variables in the GA
function store_values(fitness, approximate_grid,fr_grid, varied_params, optionsGA)%,data)

    persistent bestFitness bestApproximateGrid bestvars bestfrGrid curfitness curvars curgrid curfrgrid data_tracker

    if isempty(bestFitness)
        bestFitness = inf;
        bestApproximateGrid = [];
        bestvars = [];
        bestfrGrid = [];
        curvars = [];
        curfitness = [];
        curgrid = [];
        curfrgrid = [];
        data_tracker = {};
    end

    if(length(curfitness) == optionsGA.PopulationSize)
        curvars = [];
        curfitness = [];
        curgrid = [];
        curfrgrid = [];
    end

    % Update the best fitness and best approximate grid if current fitness is better
    if fitness < bestFitness
        bestFitness = fitness;
        bestApproximateGrid = approximate_grid;
        bestfrGrid = fr_grid;
        bestvars = varied_params;
    end
    

    curvars = [curvars; varied_params];
    curfitness = [curfitness; fitness];

    %Not sure what I was tracking here
    %data_tracker{end+1} = data;

    if isempty(curgrid)
        curgrid(:,:,end) = approximate_grid;
        curfrgrid(:,:,end) = fr_grid;
    else
        curgrid(:,:,end+1) = approximate_grid;
        curfrgrid(:,:,end+1) = fr_grid;
    end


    % Save the variables back to the main workspace
    assignin('base', 'bestFitness', bestFitness);
    assignin('base', 'bestApproximateGrid', bestApproximateGrid);
    assignin('base', 'bestvars', bestvars);
    assignin('base', 'bestfrGrid', bestfrGrid);
    assignin('base', 'curvars', curvars);
    assignin('base', 'curfitness', curfitness);
    assignin('base', 'curgrid', curgrid);
    assignin('base', 'curfrgrid', curfrgrid);
    assignin('base', 'data_tracker', data_tracker);

end