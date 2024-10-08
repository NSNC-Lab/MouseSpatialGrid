% Custom mutation function
function mutant = mutate_population(parents, optionsGA, nvars, FitnessFcn, state, thisScore, thisPopulation)

    %MUTATION


    %Was previously 0.1 and 0.01

    %Hyper Params
    %mutationRate = 0.5;
    %mutationStrength = 0.1;
    mutationRate = 0.114;  %0.2 nominal for 4D.
    mutationStrength = 0.1;

    %Grab the parents in the population.
    mutant = thisPopulation(parents,:);
    
    for i = 1:length(parents)
        cur_parent = mutant(i,:);
        for j = 1:length(cur_parent)
            if rand < mutationRate
                %mutation_out = mutant(i,j) + ((rand - 0.5) * mutationStrength);
                %Idea 1:
                mutation_out = mutant(i,j) + ((rand - 0.5) * mutationStrength);
                if mutation_out > 0
                    mutant(i,j) = mutation_out;
                end
            end
        end
    end
    

end

