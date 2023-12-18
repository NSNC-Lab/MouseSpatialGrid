function [subplot_locs,n_rows,n_cols] = detSubplotLocs(popNames)

% determines subplot locations based on cells present

% check if X cells are included
hasX = any(contains(popNames,'X'));
hasPVconverge = any(contains(popNames,'S') & ~contains(popNames,'On'));

n_rows = 3;
%n_cols = 3 + hasPVconverge + hasX;
n_cols = 4 + hasX;

topPops = {'C','R2On','R2Off','TD','S2On','S2Off','S2'};
middlePops = {'S','S1','S1On','S1Off','R1On','R1Off','ROn','ROff','SOn','SOff','X'};
bottomPops = {'On','Off'};

subplot_locs(ismember(popNames,bottomPops)) = (n_rows-1)*n_cols;
subplot_locs(ismember(popNames,topPops)) = 0;
subplot_locs(ismember(popNames,middlePops)) = (n_rows-2)*n_cols;

% all On neurons (including C) will be on the first and second columns
subplot_locs(contains(popNames,'On') | strcmp(popNames,'C')) = subplot_locs(contains(popNames,'On') | strcmp(popNames,'C')) + 2;
% PV-On should be on the first column
subplot_locs(contains(popNames,'On') & contains(popNames,'S')) = subplot_locs(contains(popNames,'On') & contains(popNames,'S')) - 1;

% Off columns should be on the third and fourth columns
subplot_locs(contains(popNames,'Off')) = subplot_locs(contains(popNames,'Off')) + n_cols;
% PV-off should be on the third
subplot_locs(contains(popNames,'Off') & contains(popNames,'S')) = subplot_locs(contains(popNames,'Off') & contains(popNames,'S')) - 1;

% for networks where there's only one PV per layer
if hasPVconverge
    subplot_locs(contains(popNames,'On') | strcmp(popNames,'C')) = subplot_locs(contains(popNames,'On') | strcmp(popNames,'C')) - 1;
    subplot_locs(contains(popNames,'Off')) = subplot_locs(contains(popNames,'Off')) - 1;
    subplot_locs(contains(popNames,'S') & ~(contains(popNames,'On') | contains(popNames,'Off'))) = ...
        subplot_locs(contains(popNames,'S') & ~(contains(popNames,'On') | contains(popNames,'Off'))) + 2;
end


% if X and TD cells are in circuit, they should be in the middle, between the on
% and off columns
subplot_locs(contains(popNames,'X') | contains(popNames,'TD')) = subplot_locs(contains(popNames,'X') | contains(popNames,'TD')) + 3;

end
