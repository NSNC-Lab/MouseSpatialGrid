function [errf_height,errf_width] = calc_errf(x,y)
% calculates the ERRF for a tuning curve. 
% INPUTS:
% x is the feature measured across e.g. sound location in degrees
% y is the mean response at each value of x e.g. mean neuronal activity at 
% each sound location

Z = trapz(x,y);
b = max(y);
a = round(Z/b);
errf_height = b;
errf_width = a;



