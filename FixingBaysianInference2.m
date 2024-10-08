warning off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% tic
% 
% 
% X_known = [1,2;
%            2,3;
%            3,4;
%            8,7;
%            2,2];
% 
% Y_known = [9;2;3;4;5];
% 
% %Get dimensional bounds
% dim_mins = min(X_known);
% dim_maxs = max(X_known);
% 
% x = linspace(dim_mins(1), dim_maxs(1), 10); % 10 points between -1 and 1 in the x-dimension
% y = linspace(dim_mins(2), dim_maxs(2), 10); % 10 points between -1 and 1 in the y-dimension
% 
% % Generate the 3D grid
% [X, Y] = ndgrid(x, y);


% X_unknown = transpose(linspace(0,0.05,2000));
% X_known = state.Population;
% Y_known = state.Score;
%X_known = [-4;-3;-2;-1;1];
%Y_known = [sin(X_known)];

%Bin size to reduce non-linearities

X_input = state.Population;
%[vals,idx] = unique(X_input);
Y_input = state.Score;

%Normalize for stability


X_known = [];
Y_known = [];

%Adjust to be the bin size at some point


max_x_width = 0.0001;
last_flag = false;

q = 1;

while(q < length(X_input))


    
    averager_x = [];
    averager_y = [];

    averager_x = [averager_x,X_input(q)];
    averager_y = [averager_y,Y_input(q)];
    
    if(q < length(X_input))
        
        if abs(X_input(q+1) - X_input(q)) < max_x_width
            
            if(q+1 == length(X_input))
                last_flag = true;
            end

            averager_x = [averager_x,X_input(q+1)];
            averager_y = [averager_y,Y_input(q+1)];
    
            start_group = true;
            counter = 1;
            while((q+counter < length(X_input)-1) && start_group == true)
                
                


                if (abs(X_input(q+1+counter) - X_input(q+counter))<max_x_width)
                    averager_x = [averager_x,X_input(q+1+counter)];
                    averager_y = [averager_y,Y_input(q+1+counter)];
                    counter = counter + 1;
                    if(q+counter == length(X_input))
                        last_flag = true;
                    end


                else
                    start_group = false;
                end
    
            end
    
            start_group = true;
    
            X_known = [X_known,mean(averager_x)];
            Y_known = [Y_known,mean(averager_y)];

            q = q+counter+1;
            %q = q + 1 + counter - 1;
        
        else
            X_known = [X_known,X_input(q)];
            Y_known = [Y_known,Y_input(q)];
            q = q+1;
        end
    

    end

    


end

if (last_flag == false)
    X_known = [X_known,X_input(q)];
    Y_known = [Y_known,Y_input(q)];
end

last_flag = false;

% n_bins = length(X_input);
% m = linspace(min(X_input),max(X_input),n_bins);
% for n =1:length(m)-1
%     bin_min = m(n);
%     bin_max = m(n+1);
% 
%     x_inter = [];
%     y_inter = [];
% 
%     for g = 1:length(X_input)
% 
%         if(X_input(g)<=bin_max && X_input(g)>bin_min)
%             x_inter = [x_inter,X_input(g)];
%             y_inter = [y_inter,Y_input(g)];
%         end
% 
%     end
% 
%     if ~anynan(mean(x_inter))
%         X_known = [X_known, mean(x_inter)];
%         Y_known = [Y_known, mean(y_inter)];
%     end
% end

X_known = transpose(X_known);
Y_known = transpose(Y_known);
% Y_known = Y_input(idx);

X_unknown = transpose(linspace(min(X_known),max(X_known),2000));
noise = 1e-2;

% X_unknown = transpose(linspace(-5,4.8,50));
% 
% X_known = transpose(-3:3);
% Y_known = sin(X_known);


% X_train = np.array([-4, -3, -2, -1, 1]).reshape(-1, 1);
% Y_train = np.sin(X_train);

% hz_search = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10];
% vert_search = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10];


% get_log_likelyhood(X_known,Y_known,0.01,0.01,noise)
% get_log_likelyhood(X_known,Y_known,0.00001,90,noise)

simplex_output = simplex_decent(0.1,0.1,X_known,Y_known,noise);
coords = simplex_output{1};
[mu,cov] = posterior_prediction(X_unknown,X_known,Y_known, coords(1), coords(2), noise);




figure;
plot(X_unknown,mu); hold on
plot(X_known,Y_known,'b.')

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
    sqdist = sum(X1.^2,2) + transpose(sum(X2.^2,2)) - 2*(X1.*transpose(X2));
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




