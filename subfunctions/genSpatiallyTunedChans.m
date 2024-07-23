function [azi,spatialCurves,chanLabels,bestLocs] = genSpatiallyTunedChans(nCells)

% generate spatially-tuned channels based on # channels
 
% +90deg is contra in mouse experiments (left), while -90 is ipsi (right)

azi = -108:108;

% fits to panniello and walker 2018, already flipped so that positive is to
% the left
load('Panniello_fits.mat','spatial_fits');
chanLabels = {spatial_fits.label};

% excitatory tuning curves
% if nCells == 1
%     chanInds = find(contains({spatial_fits.label},'center'));
%     bestLocs = 0;
%     %spatialCurves(1,:) = 0.8*gaussmf(azi,[24 0]) + 0.2; % center
% elseif nCells == 2
%     chanInds = find(contains(chanLabels,{'contra','center'}));
%     bestLocs = [90 0];
% elseif nCells == 3
%     chanInds = find(contains(chanLabels,{'contra','center','ipsi'}));
%     bestLocs = [90 0 -90];
%     %spatialCurves(1,:) = sigmf(azi,[0.1 22.5]); % contra-tuned sigmoid
%     %spatialCurves(2,:) = 0.8*gaussmf(azi,[24 0]) + 0.2; % center
%     %spatialCurves(3,:) = sigmf(azi,[-0.1 -22.5]); % ipsi-tuned sigmoid
% elseif nCells == 4
%     chanInds = find(contains(chanLabels,{'contra','45pref','center','ipsi'}));
%     bestLocs = [90 45 0 -90];
% end
% chanLabels = chanLabels(chanInds);
% 
% for i = 1:nCells
%     c_ind = find(strcmp({spatial_fits.label},chanLabels(i)));
%     spatialCurves(i,:) = spatial_fits(c_ind).curve / max(spatial_fits(c_ind).curve);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Insert spatial tuning curves from data

%Examples of implementation that runs just using gaussians

% spatialCurves(1,:) = exp(((-1/2)*(azi-90).^2)*(1/100));
% spatialCurves(2,:) = exp(((-1/2)*(azi-45).^2)*(1/100));
% spatialCurves(3,:) = exp(((-1/2)*(azi-0).^2)*(1/100));
% spatialCurves(4,:) = exp(((-1/2)*(azi+90).^2)*(1/100));

%% SOM tuning curves
contra_avg = [1	0.711835295000650	0.572164835450682	0.543664422258443];
forty5_avg = [0.769361083556803	1	0.578969909445630	0.475891141620914];
center_avg = [0.690177636848796	0.705966940305976	1	0.616779885165308];
ipsi_avg = [0.545588302488780	0.612305937959278	0.635577725387972	1];

% contra_avg = [1	0.3	0.1	0];
% forty5_avg = [0.3	1	0.25	0.2];
% center_avg = [0.3	0.3	1	0.2];
% ipsi_avg = [0.25	0.25	0.28	1];

curves = {contra_avg, forty5_avg, center_avg, ipsi_avg};
curve_data = {};
bestLocs = [90 45 0 -90];

for i=1:length(curves)
    points = [];
    curr_values = curves{i};
    for j=1:3
        rise = curr_values(j+1)-curr_values(j);
        run = bestLocs(j) - bestLocs(j+1);
        slope = rise/run;
        points_temp = (0:run-1)*slope + curr_values(j);
        points = [points points_temp];
    end
    points = [points curves{i}(end)];
    curve_data{end+1} = points;
end

spatialCurves(1,:) = curve_data{1};
spatialCurves(2,:) = curve_data{2};
spatialCurves(3,:) = curve_data{3};
spatialCurves(4,:) = curve_data{4};
azi = fliplr(linspace(-90,90, 181));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_flag = 0;

if(plot_flag == 1)

figure;
plot(azi,spatialCurves','linewidth',1);
chanNums = cellstr(string((1:nCells)'));
set(gca,'xdir','reverse');
xlabel('stimulus location'); xtickformat('degrees');
title('channel tuning curves','fontweight','normal');
xlim([-108 108]);

hold on;
locs = [90 45 0 -90];
for i = 1:4
    line([1 1]*locs(i),[0 1],'color','k','linestyle','--','displayname','')
end
legend(strcat(chanNums,{':'},chanLabels'),'location','best');


end

end