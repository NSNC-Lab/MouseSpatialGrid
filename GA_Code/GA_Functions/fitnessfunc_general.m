% Fitness function
function fitness = fitnessfunc_general(Target_grid,Target_fr_grid, varied_params, optionsGA, plot_all, all_vars,XC_vars_locs)
    
    

    %For running thing smultiple times
    %state327_18 = load('327_CH18_AbstractRun.mat');
    %varied_params = state327_18.state.bestvars{end};
    
    %First start by assigning all of the varied_params

    %Assign params by default to an area
        %Priority 1: Cross channel SOM inhibition - 3 for each hotspot
        %Priority 2: Within Channel PV wieghts - 4  total
        %Priority 3: R to C weights - 4 total
        %Priority 4: On to Ron weights - 4 total

    
    %If weights have not been set yet, intilize tem to some default value
    addpath('GA_Code\')
    cur_weights = Initilize_weights;

    %Set the weights according to 
    
    %In the base area get calculate the number of weights that you want to
    %seek (assign PV SST and hotspots and what not)
    %Check out source of noise

    %Just going to have to use all of the connections for now. Basically we
    %are running on all of the weights 

    %1. Look at source of noise
    %2. Look at optimizing the mutation
    %3. Look at wiring up weights
    %Add absolute firing rate term.
    %Add refractory to the poisson model
    %Fano factor


    %Set the weights
    %%%%%%%%%%%%%%


    counter = 0;
    for jj = 1:4
        
        if ismember(jj,XC_vars_locs)
            if jj == 1
                cur_weights(1,:) = [0,varied_params((counter*3)+1),varied_params((counter*3)+2),varied_params((counter*3)+3)];
            elseif jj == 2
                cur_weights(2,:) = [varied_params((counter*3)+1),0,varied_params((counter*3)+2),varied_params((counter*3)+3)];
            elseif jj == 3
                cur_weights(3,:) = [varied_params((counter*3)+1),varied_params((counter*3)+2),0,varied_params((counter*3)+3)];
            elseif jj == 4
                cur_weights(4,:) = [varied_params((counter*3)+1),varied_params((counter*3)+2),varied_params((counter*3)+3),0];
            end
            counter = counter + 1;
        end    
    
    end

    if all_vars(2) > 0
        cur_weights(5,:) = [varied_params(all_vars(1)+1),varied_params(all_vars(1)+2),varied_params(all_vars(1)+3),varied_params(all_vars(1)+4)];
    end

    if all_vars(3) > 0
        cur_weights(6,:) = [varied_params(sum(all_vars(1:2))+1),varied_params(sum(all_vars(1:2))+2),varied_params(sum(all_vars(1:2))+3),varied_params(sum(all_vars(1:2))+4)];
    end

    if all_vars(4) > 0
        cur_weights(7,:) = [varied_params(sum(all_vars(1:3))+1),varied_params(sum(all_vars(1:3))+2),varied_params(sum(all_vars(1:3))+3),varied_params(sum(all_vars(1:3))+4)];
    end

    GsynHandler;
    NetconHandler;

    %Run network
    %%%%%

    approximate_grid = zeros(5,4);
    fr_grid = zeros(5,4);
    
    SpikingNetwork_withOffset;
    
    %Calculate fitness
    %%%%%
    
    lambdas = [1,1,1,20,30,1];
    
    %Find hotspot rows: NOTE (XR_vars use columns see APAN poster) but
    %tuning curves occur along rows.
    hotspotrows = find(max(transpose(Target_grid(2:5,:)))>70);

    %Term 1: "Lasso Regularization"
    term1 = sum(varied_params)*lambdas(1);
    %Term 2: Normalized firing rate diff shape (clean)
    term2 = sum(abs(Target_fr_grid(1,:)./max(Target_fr_grid(1,:)) - fr_grid(1,:)./max(fr_grid(1,:))))*lambdas(2);
    if min(min(fr_grid(1,:))) == 0
        term2 = 1000;
    end
    %Term 3: Normalized firing rate diff shape (Mixed Hotspots)
    term3 = 0;
    for t3 = hotspotrows
        term3 = term3+sum(abs(Target_fr_grid(t3+1,:)./max(Target_fr_grid(t3+1,:)) - fr_grid(t3+1,:)./max(fr_grid(t3+1,:))));
    end
    if min(min(fr_grid(hotspotrows+1,:))) == 0
        term3 = 1000;
    end
    term3 = term3*lambdas(3);
    %Term 4: Performance (Clean) %10/20 Normalized
    roi_hot_idxs = find(Target_grid(1,:)>70); %Find areas of interest (Hotspots)
    term4 = sum(abs(Target_grid(1,roi_hot_idxs)-approximate_grid(1,roi_hot_idxs)))*lambdas(4);
    term4 = term4/max(Target_grid(1,:));
    %Term 5: Performance (Mixed Hotspots) %10/20 Normalized
    term5 = 0;
    for t5 = hotspotrows
        roi_hot_idxs = find(Target_grid(t5+1,:)>70); 
        term5 = term5+sum(abs(Target_grid(t5+1,roi_hot_idxs)-approximate_grid(t5+1,roi_hot_idxs)));
    end
    term5 = term5*lambdas(5);
    term5 = term5/max(max(Target_grid(2:5,:)));  %Normalize
    %Term 6: Controlling for network stability (don't let the firing rate get much higher than the ground truth)
        %Try (Don't let any point get above 30hz above the highest firing rate in the ground truth)
    term6 = max(max(Target_fr_grid))+30 - max(max(fr_grid));
    if term6 < 0
        term6 = 1000;
    else
        term6 = 0;
    end
    term6 = term6*lambdas(6);

    fitness = term1+term2+term3+term4+term5+term6;

    store_values(fitness, approximate_grid,fr_grid, varied_params,optionsGA);


    % %Start by just looking at the R to C connections
    % 
    % %R to C connections
    % varied_struct.RtoC1 = varied_params(1);
    % varied_struct.RtoC4 = varied_params(2);
    % 
    % 
    % % Assume this function updates approximate_grid
    % study_dir = fullfile(pwd,'run','4-channel-PV-inputs');
    % if exist(fullfile(study_dir, 'solve',['IC_spks_on' '.mat'])) == 2
    %     load(fullfile(study_dir, 'solve',['IC_spks_on' '.mat']),'spks','dt')
    % end
    % 
    % 
    % %Preallocate so that the appoximate grid can be populated
    % approximate_grid = zeros(5,4);
    % fr_grid = zeros(5,4);
    % 
    % SpikingNetwork_withOffset;
    % 
    % 
    % %Goal here is to match firing rates
    % fitness = abs(Target_fr_grid(1,1)-fr_grid(1,1));
    % 
    % %disp(['fitness: 1', num2str(fitness1),'  fitness: 2', num2str(fitness2),'   fitness: 3', num2str(fitness3),'   fitness: 4', num2str(fitness4)])
    % 
    % store_values(fitness, approximate_grid,fr_grid, varied_params,optionsGA);

end
