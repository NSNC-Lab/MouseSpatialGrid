%See how many hotspots there are and where to implement xc channels
XC_var_locs = find(max(Target_grid(2:5,:))>70);

%Select within channel variables to optimize
Within_Channel_PV = true;
RC_Connections = true;
On_ROn_Connections = true;


all_vars = zeros(1,4);

if ~isempty(XC_var_locs)
    all_vars(1) = length(XC_var_locs)*3;
end

if RC_Connections == true
    all_vars(2) = 4;
end

if On_ROn_Connections == true
    all_vars(3) = 4;
end

if Within_Channel_PV == true
    all_vars(4) = 4;
end

nVars = sum(all_vars);

