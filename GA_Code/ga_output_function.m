function [state,options,optchanged] = ga_output_function(options,state,flag)
persistent h1 history r
optchanged = false;
disp('Made it now')
switch flag
    case 'init'
        h1 = figure;
        ax = gca;
        ax.XLim = [0 21];
        ax.YLim = [0 21];
        l1 = min(state.Population(:,1));
        m1 = max(state.Population(:,1));
        l2 = min(state.Population(:,2));
        m2 = max(state.Population(:,2));
        r = rectangle(ax,'Position',[l1 l2 m1-l1 m2-l2]);
        history(:,:,1) = state.Population;
        assignin('base','gapopulationhistory',history);
    case 'iter'
        % Update the history every 10 generations.
        if rem(state.Generation,10) == 0
            ss = size(history,3);
            history(:,:,ss+1) = state.Population;
            assignin('base','gapopulationhistory',history);
        end
        % Find the best objective function, and stop if it is low.
        ibest = state.Best(end);
        ibest = find(state.Score == ibest,1,'last');
        bestx = state.Population(ibest,:);
        bestf = gaintobj(bestx);
        if bestf <= 0.1
            state.StopFlag = 'y';
            disp('Got below 0.1')
        end
        % Update the plot.
        figure(h1)
        l1 = min(state.Population(:,1));
        m1 = max(state.Population(:,1));
        l2 = min(state.Population(:,2));
        m2 = max(state.Population(:,2));
        r.Position = [l1 l2 m1-l1 m2-l2];
        pause(0.1)
        % Update the fraction of mutation and crossover after 25 generations.
        if state.Generation == 25
            options.CrossoverFraction = 0.8;
            optchanged = true;
        end
    case 'done'
        % Include the final population in the history.
        ss = size(history,3);
        history(:,:,ss+1) = state.Population;
        assignin('base','gapopulationhistory',history);
end
    % % Retrieve the best fitness and best approximate grid from the workspace
    % bestFitness = evalin('base', 'bestFitness');
    % bestApproximateGrid = evalin('base', 'bestApproximateGrid');
    % bestvars = evalin('base', 'bestvars');
    % 
    % %Keep track of all vars and fitness values
    % curvars = evalin('base', 'curvars');
    % curfitness = evalin('base', 'curfitness');
    % 
    % %Store the best fitness and best approximate grid after each generation
    % if ~isfield(state, 'bestFitnessHistory')
    %    gaHist.bestFitnessHistory = [];
    %    gaHist.bestApproximateGridHistory = {};
    %    gaHist.bestvarsHistory = {};
    %    gaHist.curvarsHistory = {};
    %    gaHist.curfitnessHistory = [];
    % end
    % gaHist.bestFitnessHistory(end+1) = bestFitness;
    % gaHist.bestApproximateGridHistory{end+1} = bestApproximateGrid;
    % gaHist.bestvarsHistory{end+1} = bestvars;
    % gaHist.curvarsHistory{end+1} = curvars;
    % gaHist.curfitnessHistory(end+1) = curfitness;
    % 
    % % Save the updated state to the base workspace
    % assignin('base', 'gaHist', gaHist);
    % 
    % 
    % disp('hello')
    %     optchanged = false;  % Default to false, change if options are modified
    % 
    % switch flag
    %     case 'init'
    %         % Initialization phase
    %         gaHist.bestFitnessHistory = [];
    %         gaHist.bestApproximateGridHistory = {};
    %         gaHist.bestvarsHistory = {};
    %         gaHist.curvarsHistory = {};
    %         gaHist.curfitnessHistory = [];
    % 
    %         % Save the initial state to the base workspace
    %         assignin('base', 'gaHist', gaHist);
    % 
    %     case 'iter'
    %         % Iteration phase
    %         % Retrieve the current best fitness and approximate grid
    %         bestFitness = state.Best(end);
    %         bestvars = state.Population(state.BestIndividual, :);
    % 
    %         % Custom data retrieval (make sure these variables are updated in the workspace)
    %         curvars = evalin('base', 'curvars');
    %         curfitness = evalin('base', 'curfitness');
    % 
    %         % Retrieve the stored history
    %         gaHist = evalin('base', 'gaHist');
    % 
    %         % Update the history
    %         gaHist.bestFitnessHistory(end+1) = bestFitness;
    %         gaHist.bestApproximateGridHistory{end+1} = bestApproximateGrid;  % Adjust as needed
    %         gaHist.bestvarsHistory{end+1} = bestvars;
    %         gaHist.curvarsHistory{end+1} = curvars;
    %         gaHist.curfitnessHistory(end+1) = curfitness;
    % 
    %         % Save the updated state to the base workspace
    %         assignin('base', 'gaHist', gaHist);
    % 
    %     case 'done'
    %         % Finalization phase
    %         disp('GA completed');
    % end
    % 
    % % Display a message (for debugging purposes)
    % disp('GA Output Function Called');
    % Custom Output function
%function [state, options, flag] = ga_output_function(options, state, flag)
end