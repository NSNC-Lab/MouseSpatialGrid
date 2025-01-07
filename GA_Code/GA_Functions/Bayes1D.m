tic;

warning off

input_holder = [];
input_holder_score = [];

for de = 1:length(state.curvars)
    for k = 1:length(state.curfitness{de})
        if state.curfitness{de}(k) < 200
            input_holder = [input_holder;state.curvars{de}(k,:)];
            input_holder_score = [input_holder_score;state.curfitness{de}(k)];
        end
    end
end






%X_input = input_holder;
%Y_input = input_holder_score;


%X_input = transpose(1:100);
%Y_input = sin(X_input/100);

%X_input = [input_holder(:,10)];
X_input = score(:,1);
Y_input = input_holder_score;

X_known = [];
Y_known = [];

%Adjust appropriately
max_radius = 0.0001;
last_flag = false;



%1. Look at the first thing in the list
%2. Take the distance to all of the other points in the list
%3. If something is within a certain distance then average both points and
%remove both (or all) from the list
%4. If there is no points close to the current point then remove the points
%from the list.



%Should work now. Might need some testing/troubleshooting
while(size(X_input,1)>0)
    
    collector = [];
    collector_y = [];
    qs = [];
    for q = 1:size(X_input,1)
        if norm(X_input(1,:)-X_input(q,:))<max_radius
            collector = [collector;X_input(q,:)];
            collector_y = [collector_y;Y_input(q,:)];
            qs = [qs, q];
        end
    end
    

    X_known = [X_known;mean(collector,1)];
    Y_known = [Y_known;mean(collector_y,1)];
    X_input(qs,:) = [];
    Y_input(qs,:) = [];



end



nan_X_known = [];
nan_Y_known = [];


for j = 1:length(Y_known)
    
    if ~isnan(Y_known(j))
        nan_X_known = [nan_X_known;X_known(j,:)];
        nan_Y_known = [nan_Y_known;Y_known(j,:)];
    end

end


X_known = nan_X_known;
Y_known = nan_Y_known;

% disp('hello')





% while(q < length(X_input))
% 
% 
% 
%     averager_x = [];
%     averager_y = [];
% 
%     averager_x = [averager_x,X_input(q)];
%     averager_y = [averager_y,Y_input(q)];
% 
%     if(q < length(X_input))
% 
%         if abs(X_input(q+1) - X_input(q)) < max_x_width
% 
%             if(q+1 == length(X_input))
%                 last_flag = true;
%             end
% 
%             averager_x = [averager_x,X_input(q+1)];
%             averager_y = [averager_y,Y_input(q+1)];
% 
%             start_group = true;
%             counter = 1;
%             while((q+counter < length(X_input)-1) && start_group == true)
% 
% 
% 
% 
%                 if (abs(X_input(q+1+counter) - X_input(q+counter))<max_x_width)
%                     averager_x = [averager_x,X_input(q+1+counter)];
%                     averager_y = [averager_y,Y_input(q+1+counter)];
%                     counter = counter + 1;
%                     if(q+counter == length(X_input))
%                         last_flag = true;
%                     end
% 
% 
%                 else
%                     start_group = false;
%                 end
% 
%             end
% 
%             start_group = true;
% 
%             X_known = [X_known,mean(averager_x)];
%             Y_known = [Y_known,mean(averager_y)];
% 
%             q = q+counter+1;
%             %q = q + 1 + counter - 1;
% 
%         else
%             X_known = [X_known,X_input(q)];
%             Y_known = [Y_known,Y_input(q)];
%             q = q+1;
%         end
% 
% 
%     end
% 
% 
% 
% 
% end

% if (last_flag == false)
%     X_known = [X_known,X_input(q)];
%     Y_known = [Y_known,Y_input(q)];
% end

last_flag = false;

mins = min(X_known);
maxs = max(X_known);

%X_known = transpose(X_known);
%Y_known = transpose(Y_known);


ranges = {};

% Define the ranges for each dimension

for h = 1:length(mins)
    ranges{end+1} = linspace(mins(h),maxs(h),200);
end
%ranges = {1:3, 4:6, 7:9, ... };  % Add more ranges for higher dimensions


grids = {};

% Generate the N-dimensional mesh grid
[grids{1:numel(ranges)}] = ndgrid(ranges{:});


X_unknown = [];

for x = 1:length(grids)
    values = grids{x};
    flat_vals = values(:);
    X_unknown = [X_unknown,flat_vals];
end

% Concatenate the grids along a new dimension
%X_unknown = cat(numel(ranges)+1, grids{:});



% Concatenate the grids along a new dimension
% meshgrid_ND = cat(4, X, Y, Z);



%X_unknown = transpose(linspace(min(X_known),max(X_known),2000));
%noise = 1e-2;
noise = 0.4;

simplex_output = simplex_decent(0.1,0.1,X_known,Y_known,noise);
coords = simplex_output{1};
[mu,cov] = posterior_prediction(X_unknown,X_known,Y_known, coords(1), coords(2), noise);




figure;

% plot(X_unknown,mu); hold on
% 
% plot(X_known,Y_known)

plot(grids{1},reshape(mu,size(grids{1})),'k-','LineWidth',3); hold on
plot(X_known,Y_known,'ko','MarkerSize',3)
plot(mean(X_known(1:10)),mean(Y_known(1:10)),'r.','MarkerSize',40)
plot(mean(X_known(end-9:end)),mean(Y_known(end-9:end)),'.', 'MarkerEdgeColor', [0.467 0.675 0.188], 'MarkerFaceColor', [0.467 0.675 0.188],'MarkerSize',40)


toc;

triangles = simplex_output{2};

figure;
for k = 1:length(triangles)
    plot(triangles(1,:,k),triangles(2,:,k)); hold on
end

%toc

function [mu,cov] = posterior_prediction(X_unknown,X_known,Y_known, kernal_length, kernal_vertical_variation, noise)
    
    %Usually noise would go here, but I am going to exclude temporarily
    %because eye does not accept N-dimentional arrays. 
    %Actually this shouldn't be a problem. Add back in noise if necessary
    K = kernel(X_known,X_known,kernal_length,kernal_vertical_variation)+((noise^2)*eye(length(X_known)));
    K_s = kernel(X_known,X_unknown,kernal_length,kernal_vertical_variation);
    K_ss = kernel(X_unknown,X_unknown,kernal_length,kernal_vertical_variation)+((noise)*eye(length(X_unknown)));
    
    K_inv = inv(K);
    
    mu = transpose(K_s)*K_inv*Y_known;

    cov = K_ss - transpose(K_s)*K_inv*K_s;

end

function dist = kernel(X1, X2, kernal_length,kernal_vertical_variation)
    
    %Make sure this is summing along the right axis when you go n-dim

    %Is this supposed to be matrix multiplication?

    %Last term
    dims = ndims(X2);
    order = [dims, 1:dims-1]; 
    X2_trans = permute(X2, order);
    last_term = 2*(X1*X2_trans);

    %Middle term
    summer = sum(X2.^2,2);
    dims = ndims(summer);
    order = [dims, 1:dims-1]; 
    mid_term = permute(summer, order);

    %First Term
    first_term = sum(X1.^2,2);

    sqdist = first_term + mid_term - last_term;
    dist = (kernal_vertical_variation.^2) * exp(-sqdist./(2*(kernal_length.^2)));

end



function likelyhood = get_log_likelyhood(X_known,Y_known,kernal_length,kernal_vertical_variation,noise)
    K = kernel(X_known,X_known,kernal_length,kernal_vertical_variation)+((noise^2)*eye(length(X_known)));
    %Transpose might be backwards
    %Might also want to pass in values like det(K) and inv(K). Expensive
    %calculates can be down before hand instead of every iteration.
    
    %Account for logorithm dicontinuity
    if det(K) < realmin
        surrogate = log(realmin);

    else
        surrogate = log(det(K));
    end

    likelyhood = (0.5*surrogate) + (0.5*transpose(Y_known)*inv(K)*Y_known) + 0.5*(length(X_known)*log(2*pi));


    %Not sure how okay this is
    if likelyhood < 0
        likelyhood = inf;
    end

    %disp(likelyhood)

    
end

function pair_and_plottables = simplex_decent(l_start,v_start,X_known,Y_known,noise)
    
    
    %Hyper Params
    l = l_start;
    v = v_start;
    simplex_mag = 0.0001;
    curgen = 1;

    Max_Gens = 1000;

    triangle_keeper = zeros(2,6,Max_Gens);
    
    saved_likelyhoods = [];

    gains = [1e-2,1e-1,1e0,1e1,1e2];

    % flag = false;
    % ingen_max = 1000;
    % ingen = 0;

    % simplex_mag_up = l;
    % simplex_mag_down = l;

    %Do initial check
    likelyhood_init = get_log_likelyhood(X_known,Y_known,l_start,v_start,noise);
    saved_likelyhoods = [saved_likelyhoods,likelyhood_init];

    while(curgen<Max_Gens)
        
        likelyhoods = [];
        ls = [];
        vs = [];
        %Just for plotting purposes track x.
        xs = [];

        for x = gains
            
            %Triangle
            %l_test = [l, l+simplex_mag*0.866*x,l-simplex_mag*0.866*x];
            %v_test = [v+simplex_mag*x,v-simplex_mag*0.5*x,v-simplex_mag*0.5*x];
        
            %Hexegon
            l_test = [l+simplex_mag*0.5*x,l+simplex_mag*x,l+simplex_mag*0.5*x,l-simplex_mag*0.5*x,l-simplex_mag*x,l-simplex_mag*0.5*x];
            v_test = [v+simplex_mag*x*0.866,v,v-simplex_mag*0.866*x,v-simplex_mag*0.866*x,v,v+simplex_mag*x*0.866];

            for z = 1:length(l_test)
                if l_test(z) < 0.000001
                    %Might need to adjust this
                    l_test(z) = 0.000001;
                end

                if v_test(z) > 1
                    v_test(z) = 1;
                end
        
                if v_test(z) < 0
                    v_test(z) = 0.0001;
                end
               
                likelyhood = get_log_likelyhood(X_known,Y_known,l_test(z),v_test(z),noise);
                likelyhoods = [likelyhoods,likelyhood];
                ls = [ls,l_test(z)];
                vs = [vs,v_test(z)];
                xs = [xs,x];


            end
            
            % likelyhood1 = get_log_likelyhood(X_known,Y_known,l_test(1),v_test(1),noise);
            % likelyhood2 = get_log_likelyhood(X_known,Y_known,l_test(2),v_test(2),noise);
            % likelyhood3 = get_log_likelyhood(X_known,Y_known,l_test(3),v_test(3),noise);
            
            
            
            
        
            %ingen = ingen + 1;
        
          
        end


        [output_likelyhood,index_in] = min(likelyhoods);

        %if output_likelyhood < saved_likelyhoods(end)
            saved_likelyhoods = [saved_likelyhoods,output_likelyhood];
            % ingen = 0;
            % flag = true;
            l = ls(index_in);
            v = vs(index_in);
            x = xs(index_in);
        %end

    % triangle_keeper(:,:,curgen) = [l+simplex_mag*0.5*x,l+simplex_mag*x,l+simplex_mag*0.5*x,l-simplex_mag*0.5*x,l-simplex_mag*x,l-simplex_mag*0.5*x];
    %           v+simplex_mag,v-simplex_mag*0.5,v-simplex_mag*0.5];
        
    %Hexegon
    triangle_keeper(:,:,curgen) = [l+simplex_mag*0.5*x,l+simplex_mag*x,l+simplex_mag*0.5*x,l-simplex_mag*0.5*x,l-simplex_mag*x,l-simplex_mag*0.5*x;
            v+simplex_mag*x*0.866,v,v-simplex_mag*0.866*x,v-simplex_mag*0.866*x,v,v+simplex_mag*x*0.866];

    curgen = curgen + 1;
            

    end

    optimal_pair = [l,v];

    pair_and_plottables = {optimal_pair, triangle_keeper};

    % flag = false;
   

    

end




