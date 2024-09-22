% Custom mutation function
function mutant = mutate_population(parents, optionsGA, nvars, FitnessFcn, state, thisScore, thisPopulation)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %Relaxing ROI
    % if(state.Generation < 4)  %We need atleast 2 recordings in the state struct. There is a bug with the first generations and tracking the fitness correclty there.
    %     mutationStrength = 0.12;
    %     mutationRate = 0.2;
    % 
    %     assignin('base', 'mutationStrength', mutationStrength);
    %     assignin('base', 'mutationRate', mutationRate);
    % else
    %     mutationStrength = evalin('base','mutationStrength');
    %     mutationRate = evalin('base','mutationRate');
    % 
    %     previous_mean = mean(state.curfitness{end - 1});
    % 
    %     % if previous_mean - 1 > mean(state.curfitness{end}) %Most stringent, must go first
    %     %     mutationStrength = mutationStrength/2;
    %     %     mutationRate = mutationRate/2;
    %     % elseif previous_mean > mean(state.curfitness{end})
    %     %     mutationStrength = mutationStrength/1.3;
    %     %     mutationRate = mutationRate/1.3;
    %     if previous_mean > mean(state.curfitness{end})
    %         mutationStrength = mutationStrength/1.5;
    %         mutationRate = mutationRate/1.05;
    %     else
    %         mutationStrength = mutationStrength*1.5;
    %         mutationRate = mutationRate*1.05;
    %     end
    % 
    %     assignin('base', 'mutationStrength', mutationStrength);
    %     assignin('base', 'mutationRate', mutationRate);
    % 
    % end
    % 
    % %Should add some sort of benefit
    % %Think about it like there are onyl so much benfit ot give. Once the
    % %algorithm has gotten all of the beneift for the benefit supply it will
    % %increase its mutation strength and mutation rate much less rapidly.
    % 
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%
    %Original
    %mutationStrength = 0.2;  
    mutationStrength = 0.2;
    mutationRate = 0.2;

    %disp(mutationRate);
    %disp(mutationStrength);
    
    %Currently mutates every parameter of a given creature.
    %Instead we should mutate along each axis, given by the importance of
    %the mutation. Each gene should have a chance to mutate instead of each
    %creature.


    mutant = thisPopulation(parents,:);
    for i = 1:size(mutant, 1)
        if rand < mutationRate
            mutant(i,:) = mutant(i,:) + (rand(1, nvars) - 0.5) * mutationStrength;
        end
    end
    %disp('Made it here')
end
