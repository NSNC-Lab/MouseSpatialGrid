% update netcons
NetconHandler;

%Number of nodes (x data)
%heirarchy of nodes (y data)
channels = [1:4];

numChannels = length(channels);
numPopulations = length(s.populations);
numConnections = length(s.connections);

% population without C
le = numPopulations-1;

% Total # of nodes, channels*(population without C) plus 1 for output C
totalNodes = (numChannels*le)+1;

% Preallocate for efficiency
Nodesx = zeros(1,totalNodes);
Nodesy = zeros(1,totalNodes);

% node counter as index
nodeCounter = 0;

%% Mapping of population names to x and y offsets
popNameOffsets = containers.Map(...
    {'On', 'Off', 'ROn', 'ROff', 'SOnOff', 'TD', 'X'}, ...
    {[0, 0], [0.5, 0], [0, 2], [0.5, 2], [0.25, 1], [0.75, 3], [-0.25, 1]});

for b = 1:numChannels
    for a = 1:numPopulations
        popName = s.populations(a).name;

        if isKey(popNameOffsets, popName)
            offsets = popNameOffsets(popName);
            nodeCounter = nodeCounter + 1;
            Nodesx(nodeCounter) = channels(b) + offsets(1);
            Nodesy(nodeCounter) = offsets(2);
        end

    end
end

%% single offset implementation for C
nodeCounter = nodeCounter + 1;
Nodesx(nodeCounter) = (numChannels+1)/2;
Nodesy(nodeCounter) = 3.5;


%% Determining the start and end node #'s
%Should be able to somehow automatically populate this with s.connections
sources = [];
targets = [];

inhibs_PEs = [];
inhibs_PEt = {};

inhibs_XRs = [];
inhibs_XRt = [];

TD_SOMs = [];
TD_SOMt = [];
[XR_m,XR_n] = size(netcons.XRnetcon);


%Two sanity checks
%If only connect contra make sure c looks like contra
%All connections rerun on normal tuning curves

%% create graphable arrays of start and end node locations for each connection
cons_possible = {'On->On', 'Off->Off', 'On->ROn', ...
                 'On->SOnOff', 'SOnOff->Ron', 'SOnOff->ROff', ...
                 'Off->ROff', 'Off->SOnOff', 'ROn->ROn', ...
                 'ROn->X', 'X->ROn', 'TD->ROn', ...
                 'TD->ROff', 'TD->X', 'ROn->C'};

inhibitory_cons = {'SOnOff', 'X', 'TD'};

% mapping for node numberings
nodeMap = containers.Map({'On', 'Off', 'ROn', 'ROff', 'SOnOff', 'TD', 'X', 'C'}, ...
                           {1, 2, 3, 4, 5, 6, 7, totalNodes});

for c = 1:numChannels
    for x = 1:numConnections
        % grab direction string
        direction = s.connections(x).direction;

        % using current netcon should make this dynamic
        netcon_index = find(strcmp('netcon', s.connections(x).parameters)) + 1;
        current_netcon = s.connections(x).parameters(netcon_index);
        current_netcon = current_netcon{1};
        [m,n] = size(current_netcon);
        
        % get start node and end node
        str = extractBefore(direction,"-");
        stra = extractAfter(direction,">");

        % Check if the connection is inhibitory
        if ismember(str, inhibitory_cons)
            for p = 1:n
                if current_netcon(c, p) == 1
                    s_inhib = (c-1) * le + nodeMap(str);
                    t_inhib = (p-1) * le + nodeMap(stra);
                    if strcmp(str, 'X')
                        inhibs_XRt = [inhibs_XRt t_inhib];
                        inhibs_XRs = [inhibs_XRs s_inhib];
                    elseif strcmp(str, 'TD') && strcmp(stra, 'X')
                        TD_SOMs = [TD_SOMs s_inhib];
                        TD_SOMt = [TD_SOMt t_inhib];
                    elseif strcmp(str, 'SOnOff')
                        if c==p % within channel, flag 1
                            inhibs_PEt{end+1} = {t_inhib, 1};
                        else    % cross channel, flag 2
                            inhibs_PEt{end+1} = {t_inhib, 2};
                        end
                        inhibs_PEs = [inhibs_PEs s_inhib];
                    end
                    targets = [targets t_inhib];
                    sources = [sources s_inhib];
                end
            end
        else
            % Handle non-inhibitory connections
            for p = 1:n
                if current_netcon(c, p) == 1
                    % assign source node a number
                    source = (c-1) * le + nodeMap(str);

                    % if target node is C, label as last node
                    if strcmp(stra, 'C')
                        tar = nodeMap(stra);
                    else
                        % assign target node a number
                        tar = (p-1) * le + nodeMap(stra);
                    end
                    % add assigned node # to source and target arrays
                    sources = [sources source];
                    targets = [targets tar];
                end
            end
        end
    end
end


%% Graphing Segment

%Be careful! Do not override the s struct
%s1 = [1 1 2 2 5 5 3 6 6 3];
%t = [3 5 4 5 3 4 8 3 4 7];
weights = ones(1,length(targets))*1;

figure;

G = digraph(sources,targets,weights);
%G = digraph(s,t)

%IMPORTANT:::: Look at
%https://www.mathworks.com/help/matlab/ref/matlab.graphics.chart.primitive.graphplot.highlight.html
%for graph colors

%% graph and change X and PV attributes
p = plot(G,'k','Xdata',Nodesx,'Ydata',Nodesy); hold on

%% highlight P-E connections
highlight(p, inhibs_PEs,'NodeColor','red')
for node=1:length(inhibs_PEt)
    if inhibs_PEt{node}{2} == 1 % within channel
        highlight(p, inhibs_PEs(node), inhibs_PEt{node}{1},'EdgeColor','blue')
    else                        % cross channel
        highlight(p, inhibs_PEs(node), inhibs_PEt{node}{1},'EdgeColor','green')
    end
end

%% highlight X-R connections
highlight(p, inhibs_XRs,'NodeColor','red')
highlight(p, inhibs_XRs, inhibs_XRt,'EdgeColor','red')
%highlight(p, inhibs_XRs, inhibs_XRt,'LineStyle','--')

%% highlight TD-X connections
highlight(p, TD_SOMs,'NodeColor','red')
highlight(p, TD_SOMs, TD_SOMt,'EdgeColor','yellow')


%% synaptic strength modifier
%Now we have to go through and "highlight" which means to change the size
%of the line widths that connect neurons.
syns = [];

%Might have to do this for inhibs_2/s2
for u = 1:length(sources)
    start_type = mod(sources(u),7);
    end_type = mod(targets(u),7);

    if start_type == 0
        start_type = 7;
    end
    if end_type == 0
        end_type = 7;
    end

    tstring = s.populations(start_type).name + "->" + s.populations(end_type).name;
    for j = 1:length(s.connections)
        if strcmp(tstring,s.connections(j).direction)
            target = j;
        end
    end

    gsyn_m = s.connections(target).parameters(2);

    if gsyn_m{1} > 0.2
        end_gsyn = 0.5;
    else
        end_gsyn = (gsyn_m{1}^1.5)*750;
    end

    syns = [syns gsyn_m];

    highlight(p, sources(u), targets(u), 'LineWidth',end_gsyn,'ArrowSize',end_gsyn*5)


end


%% GPT HSV test
% syns = [];
% 
% for u = 1:length(sources)
%     start_type = mod(sources(u),7);
%     end_type = mod(targets(u),7);
% 
%     if start_type == 0
%         start_type = 7;
%     end
%     if end_type == 0
%         end_type = 7;
%     end
% 
%     tstring = s.populations(start_type).name + "->" + s.populations(end_type).name;
%     for j = 1:length(s.connections)
%         if strcmp(tstring,s.connections(j).direction)
%             target = j;
%         end
%     end
% 
%     gsyn_m = s.connections(target).parameters(2);
%     syns = [syns gsyn_m{1}];
% 
% end
% 
% % Step 2: Normalize the gsyn values
% min_gsyn = min(syns);
% max_gsyn = max(syns);
% norm_gsyns = (syns - min_gsyn) / (max_gsyn - min_gsyn);  % Normalize to [0, 1]
% 
% % Step 3: Assign hues based on normalized gsyn values
% for u = 1:length(sources)
%     start_type = mod(sources(u), 7);
%     end_type = mod(targets(u), 7);
% 
%     if start_type == 0
%         start_type = 7;
%     end
%     if end_type == 0
%         end_type = 7;
%     end
% 
%     tstring = s.populations(start_type).name + "->" + s.populations(end_type).name;
%     for j = 1:length(s.connections)
%         if strcmp(tstring, s.connections(j).direction)
%             target = j;
%         end
%     end
% 
%     gsyn_m = s.connections(target).parameters(2);
%     gsyn_value = gsyn_m{1};
% 
%     % Find the normalized value for this gsyn
%     norm_value = (gsyn_value - min_gsyn) / (max_gsyn - min_gsyn);
% 
%     % Set hue based on normalized gsyn value
%     hue = norm_value;  % Hue varies from 0 to 1
%     saturation = 1;    % Full saturation
%     value = 1;         % Full brightness
%     color = hsv2rgb([hue, saturation, value]);
% 
%     % Highlight the line with the calculated color
%     highlight(p, sources(u), targets(u), 'EdgeColor', color, 'ArrowSize', 5);
% end


%What are the things that go into the magnitude of the connection?
%Gsyn
%Netcon
    %PEnetcon
    %RCnetcon

%annotation('arrow', [.56 .58], [.221 .223],'color',[0 0 1]);

%ToDos
%Make all of the connections dynamic    %Still needs some work
      %Should make it so that if a new node  is added to the system or
      %taken away it still works.
      %Should work with different cortical neuron setups.

%Extend to 4 channels                   %Done
%Change colors                          %Mostly done
%Change line/arrow weights              %Mostly done
%Put in arrows                          %Done
%Find a way to handle x inhibition without it being too crowded   %Done

%Dales law implement   %Mostly Done

%Things to touch on in today's meeting
%1. Netcon is not implemented yet for line width multiplier
%2. Labels on lines
%3. Further flexibility

%Todos (After meeting)
%Make TDs red
%Replicate Aim network model