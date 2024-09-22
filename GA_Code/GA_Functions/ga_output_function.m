% Custom Output function
function [state, optionsGA, optchanged] = ga_output_function(optionsGA, state, flag)
    optchanged = false;

    % Retrieve the best fitness and best approximate grid from the workspace
    bestFitness = evalin('base', 'bestFitness');
    bestApproximateGrid = evalin('base', 'bestApproximateGrid');
    bestvars = evalin('base', 'bestvars');
    bestfrGrid = evalin('base', 'bestfrGrid');
    curfitness = evalin('base','curfitness');
    curvars = evalin('base','curvars');
    curgrid = evalin('base','curgrid');
    curfrgrid = evalin('base','curfrgrid');

    %Store the best fitness and best approximate grid after each generation
    if ~isfield(state, 'bestFitnessHistory')
       state.bestFitnessHistory = [];
       state.bestApproximateGridHistory = {};
       state.bestvars = {};
    
       state.frGridHistory = {};
       state.curfitness = {};
       state.curvars = {};
       state.curgrid = {};
       state.curfrgrid = {};
       state.convergence = [];
    end
    
    if(state.Generation >0)
        state.bestFitnessHistory(end+1) = bestFitness;
        state.bestApproximateGridHistory{end+1} = bestApproximateGrid;
        state.bestvars{end+1} = bestvars;

        state.frGridHistory{end+1} = bestfrGrid;
        state.curfitness{end+1} = curfitness;
        state.curvars{end+1} = curvars;
        state.curgrid{end+1} = curgrid;
        state.curfrgrid{end+1} = curfrgrid;

        %Implement convergence quantification here.
    
        %Find centroid
        centroid = mean(state.Population); 
        %Calculate distance and sum each point
        dist_mat = sum(sum(abs(state.Population-centroid)));
        %Track
        state.convergence(end+1) = dist_mat;

        %FixingBaysianInferenceND;
        %state.Population = [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];

    end

    % Save the updated state to the base workspace
    assignin('base', 'state', state);
end