nodeNames = {'In','E','PV','SOM','VIP'};

% 5 directed edges: s -> t
s = [1 2 2 4 5]';
t = [2 3 4 2 2]';

% weights from your optimizer (length = 5)
w = [0.8 -0.4 0.2 -0.7 0.1]';

% manual positions (stable + clean)
xy = [0 0; 1 0; 2 0.7; 2 -0.7; 3 0];   % N x 2

plot_small_network(nodeNames, s, t, w, xy);

function h = plot_small_network(nodeNames, s, t, w, xy)
% nodeNames : 1xN cellstr
% s,t       : Ex1 integer node indices (1..N)
% w         : Ex1 weights (line widths/colors derived from these)
% xy        : Nx2 node positions [x y] (recommended for consistent layout)

    arguments
        nodeNames (1,:) cell
        s (:,1) double
        t (:,1) double
        w (:,1) double
        xy (:,2) double
    end

    N = numel(nodeNames);

    % Directed graph (use graph(...) if undirected)
    G = digraph(s, t, w, nodeNames);

    % Map |w| -> line width (tune these)
    lwMin = 0.5; lwMax = 6;
    a = abs(w);
    a0 = min(a); a1 = max(a);
    lw = lwMin + (a - a0) ./ (a1 - a0 + eps) * (lwMax - lwMin);

    figure('Color','w');
    h = plot(G, ...
        'XData', xy(:,1), 'YData', xy(:,2), ...
        'NodeLabel', nodeNames, ...
        'ArrowSize', 14, ...
        'MarkerSize', 10, ...
        'LineWidth', lw);

    % Edge color mapped by weight value (diverging look)
    h.EdgeCData  = w;
    h.EdgeColor  = 'flat';
    colormap(parula);                 % swap if you want a diverging map
    c = colorbar; c.Label.String = 'weight';

    % Cosmetics
    h.NodeFontSize = 12;
    h.NodeFontWeight = 'bold';
    axis off;
end

