% Custom population creation function
function population = create_population(GenomeLength, ~, optionsGA)
    
    %%%%  OPTION 1 %%%%%    (Random initilization)
    n = optionsGA.PopulationSize;
    %population = rand(n, GenomeLength)*0.15;
    
    %More optimal start

    %Shift to the right
    %population(:,1) = population(:,1) + 0.15;

    %a = fieldnames(optionsGA);

    %More optimal start
    population = rand(n, GenomeLength)*0.05;
    % population(:,3) = population(:,3)*2;
    % population(:,2) = population(:,2)*6;
    

    % for k = 1:length(a)
    %     if strcmp(a{k},'InitialPopulationMatrix')
    %         disp('hello')
    %         population = optionsGA.InitialPopulationMatrix
    %     end
    % end

    %%%%  OPTION 2 %%%%%    (Evenly spaced starting)
    % %Try an evenly spaced grid
    % N = 50; % desired total number of points
    % range = [0, 0.15]; % range for both x and y directions
    % 
    % % Define the number of points and range
    % % N: desired total number of points
    % % range: 2-element vector [min, max] for both x and y directions
    % 
    % % Calculate the spacing needed to fit N points in the given range
    % area = (range(2) - range(1))^2;
    % d = sqrt(area / (N * sqrt(3)/2));
    % 
    % % Determine number of points in each direction
    % nx = ceil((range(2) - range(1)) / d);
    % ny = ceil((range(2) - range(1)) / (d * sqrt(3)/2));
    % 
    % % Create hexagonal grid
    % X = [];
    % Y = [];
    % for j = 0:ny
    %     for i = 0:nx
    %         x = range(1) + i * d;
    %         y = range(1) + j * d * sqrt(3)/2;
    %         if mod(j, 2) == 1
    %             x = x + d/2;
    %         end
    %         if x <= range(2) && y <= range(2)
    %             X = [X, x];
    %             Y = [Y, y];
    %         end
    %     end
    % end
    % 
    % % Trim to desired number of points if necessary
    % if length(X) > N
    %     idx = randperm(length(X), N);
    %     X = X(idx);
    %     Y = Y(idx);
    % end
    % 
    % population = transpose([X;Y]);
    % 
    % %disp('hello')
end