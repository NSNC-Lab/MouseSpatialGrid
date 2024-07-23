% Custom Output function
function [state, optionsGA, optchanged] = ga_output_function(optionsGA, state, flag)
    optchanged = false;

    % Retrieve the best fitness and best approximate grid from the workspace
    bestFitness = evalin('base', 'bestFitness');
    bestApproximateGrid = evalin('base', 'bestApproximateGrid');
    bestvars = evalin('base', 'bestvars');
    curfitness = evalin('base','curfitness');
    curvars = evalin('base','curvars');
    curgrid = evalin('base','curgrid');

    %Store the best fitness and best approximate grid after each generation
    if ~isfield(state, 'bestFitnessHistory')
       state.bestFitnessHistory = [];
       state.bestApproximateGridHistory = {};
       state.bestvars = {};
       state.curfitness = {};
       state.curvars = {};
       state.curgrid = {};
    end
    
    if(state.Generation >1)
        state.bestFitnessHistory(end+1) = bestFitness;
        state.bestApproximateGridHistory{end+1} = bestApproximateGrid;
        state.bestvars{end+1} = bestvars;
        state.curfitness{end+1} = curfitness;
        state.curvars{end+1} = curvars;
        state.curgrid{end+1} = curgrid;
    end


    % Save the updated state to the base workspace
    assignin('base', 'state', state);
end