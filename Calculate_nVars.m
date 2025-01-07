%See how many hotspots there are and where to implement xc channels
XC_var_locs = max(transpose(Target_grid(2:5,:)))>70;

%Select within channel variables to optimize
RC_Connections = false;
On_ROn_Connections = false;
Within_Channel_PV = false;


all_vars = zeros(1,4);

if ~isempty(X_var_locs)
    all_vars(1) = length(XC_var_locs)*3;
end

if RC_Connections == true
    all_vars(2) = 4;
end

if On_ROn_Connections == true
    all_vars(3) = 4;
end

if RC_Connections == true
    Within_Channel_PV(4) = 4;
end

nVars = sum(all_vars);

