test_netcons = containers.Map();

for x=1:length(s.connections)
    direction = s.connections(x).direction;
    netcon_index = find(strcmp('netcon', s.connections(x).parameters)) + 1;
    current_netcon = s.connections(x).parameters{netcon_index};
    test_netcons(direction) = current_netcon;
end

disp(test_netcons)
