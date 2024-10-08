

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


X_input = X_input(1:3);
Y_input = Y_input(1:3);


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
noise = 1e-8;

% X_unknown = transpose(linspace(-5,4.8,50));
% 
% X_known = transpose(-3:3);
% Y_known = sin(X_known);


% X_train = np.array([-4, -3, -2, -1, 1]).reshape(-1, 1);
% Y_train = np.sin(X_train);

% hz_search = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10];
% vert_search = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10];


simplex_output = simplex_decent(0.000001,0.000001,X_known,Y_known,noise);
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
    simplex_mag = 0.01;
    curgen = 1;

    Max_Gens = 1000;

    triangle_keeper = zeros(2,3,Max_Gens);
    
    saved_likelyhoods = [];



    flag = false;
    ingen_max = 1000;
    ingen = 0;

    simplex_mag_up = l;
    simplex_mag_down = l;

    %Do initial check
    likelyhood_init = get_log_likelyhood(X_known,Y_known,l_start,v_start,noise);
    saved_likelyhoods = [saved_likelyhoods,likelyhood_init];

    while(curgen<Max_Gens)
        

        



        %Pick up tommorow
        %Do scount hikes around the space
            %Increase it and decrease it iteratively until you find a
            %better place to go.
                %If not better place then stop.
        %while(flag == false && ingen < ingen_max)
            %Start exploring
            %Do inital check
            % if ingen == 0

                
                l_backup = l;
                v_backup = v;
                
                

                l = [l, l+simplex_mag*0.866,l-simplex_mag*0.866];
                v = [v+simplex_mag,v-simplex_mag*0.5,v-simplex_mag*0.5];

                for z = 1:length(l)
                    if l(z) < 0
                        %Might need to adjust this
                        l(z) = 0.0001;
                    end

                    if v(z) < 0
                        v(z) = 0.0001;
                    end
                end

                likelyhood1 = get_log_likelyhood(X_known,Y_known,l(1),v(1),noise);
                likelyhood2 = get_log_likelyhood(X_known,Y_known,l(2),v(2),noise);
                likelyhood3 = get_log_likelyhood(X_known,Y_known,l(3),v(3),noise);
                [output_likelyhood,index_in] = min([likelyhood1,likelyhood2,likelyhood3]);

                ingen = ingen + 1;

                if output_likelyhood < saved_likelyhoods(end)
                    saved_likelyhoods = [saved_likelyhoods,output_likelyhood];
                    ingen = 0;
                    flag = true;
                    l = l(index_in);
                    v = v(index_in);
                end

                

            
            %else
                % %Scale up
                % 
                % simplex_mag_up = simplex_mag_up * 1.05;
                % 
                % l_up = [l_backup, l_backup+simplex_mag_up*0.866,l_backup-simplex_mag_up*0.866];
                % v_up = [v_backup+simplex_mag_up,v_backup-simplex_mag_up*0.5,v_backup-simplex_mag_up*0.5];
                % 
                % for z = 1:length(l)
                %     if l_up(z) < 0
                %         l_up(z) = 0.0001;
                %     end
                % 
                %     if v_up(z) < 0
                %         v_up(z) = 0.0001;
                %     end
                % end
                % 
                % likelyhood1u = get_log_likelyhood(X_known,Y_known,l_up(1),v_up(1),noise);
                % likelyhood2u = get_log_likelyhood(X_known,Y_known,l_up(2),v_up(2),noise);
                % likelyhood3u = get_log_likelyhood(X_known,Y_known,l_up(3),v_up(3),noise);
                % 
                % simplex_mag_down = simplex_mag_down * 0.95;
                % 
                % l_down = [l_backup, l_backup+simplex_mag_down*0.866,l_backup-simplex_mag_down*0.866];
                % v_down = [v_backup+simplex_mag_down,v_backup-simplex_mag_down*0.5,v_backup-simplex_mag_down*0.5];
                % 
                % for z = 1:length(l)
                %     if l_down(z) < 0
                %         l_down(z) = 0.0001;
                %     end
                % 
                %     if v_down(z) < 0
                %         v_down(z) = 0.0001;
                %     end
                % end
                % 
                % %Scale down
                % likelyhood1d = get_log_likelyhood(X_known,Y_known,l_down(1),v_down(1),noise);
                % likelyhood2d = get_log_likelyhood(X_known,Y_known,l_down(2),v_down(2),noise);
                % likelyhood3d = get_log_likelyhood(X_known,Y_known,l_down(3),v_down(3),noise);
                % 
                % [output_likelyhood,index_in] = min([likelyhood1u,likelyhood2u,likelyhood3u,likelyhood1d,likelyhood2d,likelyhood3d]);
                % 
                % 
                % ingen = ingen + 1;
                % 
                % if output_likelyhood < saved_likelyhoods(end)
                %     saved_likelyhoods = [saved_likelyhoods,output_likelyhood];
                %     ingen = 0;
                %     flag = true;
                %     if index_in <= 3
                %         simplex_mag = simplex_mag_up;
                %         simplex_mag_down = simplex_mag_up;
                %         l = l_up(index_in);
                %         v = v_up(index_in);
                % 
                %     else
                %         simplex_mag = simplex_mag_down;
                %         simplex_mag_up = simplex_mag_down;
                %         l = l_down(index_in-3);
                %         v = v_down(index_in-3);
                %     end
                % 
                % end

                

           % end



            

        end

        flag = false;


        %Get points around current point
        
        
       
        triangle_keeper(:,:,curgen) = [l, l+simplex_mag*0.866,l-simplex_mag*0.866;
                      v+simplex_mag,v-simplex_mag*0.5,v-simplex_mag*0.5];
        
        % l = [l, l+simplex_mag*0.866,l-simplex_mag*0.866];
        % v = [v+simplex_mag,v-simplex_mag*0.5,v-simplex_mag*0.5];
        % 
        % %Search point 1
        % likelyhood1 = get_log_likelyhood(X_known,Y_known,l(1),v(1),noise);
        % likelyhood2 = get_log_likelyhood(X_known,Y_known,l(2),v(2),noise);
        % likelyhood3 = get_log_likelyhood(X_known,Y_known,l(3),v(3),noise);
        % 
        % 
        % 
        % [likelymin,index] = min([likelyhood1,likelyhood2,likelyhood3]);

        % if curgen>2
        %     if likelymin > saved_likelyhoods(end-1)
        %         simplex_mag = simplex_mag/1.5;
        %     end
        % end
        
       

        curgen = curgen + 1;

        % saved_likelyhoods = [saved_likelyhoods,likelymin];
    
    % end
    
    
    optimal_pair = [l,v];

    pair_and_plottables = {optimal_pair, triangle_keeper};
    

end


% function best_mu = grid_serach_params(sigma_hz,sigma_vert,X_unknown,X_known,Y_known,noise)
% 
% 
%     best_parameter_pair = inf;
%     best_i = 0;
%     best_j = 0;
% 
%     for i = 1:length(sigma_hz)
%         for j = 1:length(sigma_vert)
% 
% 
% 
%             [mu,cov] = posterior_prediction(X_unknown,X_known,Y_known, sigma_hz(i), sigma_vert(j), noise);
% 
% 
%             residuals = [];
% 
%             for k = 1:size(X_known)
% 
%                 % Find the closest x-coordinate on the real line
%                 [~, closest_idx] = min(abs(mu(:, 1) - Y_known(k)));
% 
%                 % Calculate the vertical (y) difference
%                 vertical_residual = abs(mu(closest_idx) - Y_known(k));
% 
%                 % Store the residual
%                 residuals = [residuals , vertical_residual];
%             end
% 
%             if sum(residuals)<best_parameter_pair
%                 best_parameter_pair = sum(residuals);
%                 best_i = i;
%                 best_j = j;
%             end
% 
% 
%             figure;
%             plot(X_unknown,mu); hold on
%             plot(X_known,Y_known,'b.')
% 
% 
%         end
%     end
% 
% 
%     [best_mu,~] = posterior_prediction(X_unknown,X_known,Y_known, sigma_hz(best_i), sigma_vert(best_j), noise);
% 
% 
% 
% 
% end

