function [c,ceq] = normFiringRate(w)

% sum of weights must not exceed # channels
c(1) = sum(w(w >= 0 & w <= 1)) - sum(w > 0);

% weights must not be negative
ceq(1) = length(w) - sum(w == 0) - sum(w > 0);

end