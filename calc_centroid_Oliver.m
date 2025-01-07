function centroid = calc_centroid_Oliver(x,y)
% Work out the centroid for each unit (see Middlebrooks & Bremen 2013). The
% centroid was computed from an RAF by finding the peak range of one or
% more contiguous sound-source locations that elicited spike rates >=75% of
% a neuron’s maximum rate plus the two locations on either side of that
% range. All the locations within the peak range were treated as vectors
% weighted by their corresponding spike rates. A vector sum was formed, and
% the direction of the resultant vector was taken as the centroid.
% Also see this website https://urldefense.com/v3/__http://www.1728.org/vectutor.htm__;!!DZ3fjg!_lZLZrdVHoTKgDfZZNHO6tzdlLEuDUJRwwKVr6tQ6TkM3LnhmpCgwqe9qFh090JHJHEo4ML9RnzauFl23NqJFvhbNUUTTw$ 

% INPUTS
% x = sound location in degrees
% y = mean response to each value of x
% plot_fig = true/false if you would like to plot the outcome or not


%% Put degree values on scale of 0-180 degrees
% with zero at 90 degrees right and a counterclockwise direction

% ****** you may need to edit this. I had my x values so that 0 degrees was
% positioned in the middle. the left hemifield was negative numbers so 90
% degrees to the left was -90 and right hemifield was positive, so 90
% degrees to the right was +90... ******************************


% orig_x = x;
% x(orig_x==0) = 90;
% x(orig_x<0) = abs(orig_x(orig_x<0))+90;
% x(orig_x>0) = abs(orig_x(orig_x>0)-90);


%% Convert to radians
xRad = zeros(size(x));
for jj = 1:length(x)
    xRad(jj) = x(jj)/180*pi;
end

%% Find data points/locations within 75% of max firing rate - these
% will contribute to the centroid
peakRange = (y>(max(y)*0.75));
y = y-min(y); % normalize y to between 0 and 1
y = y/max(y)*100; % magnitudes, see this website https://urldefense.com/v3/__http://www.1728.org/vectutor.htm__;!!DZ3fjg!_lZLZrdVHoTKgDfZZNHO6tzdlLEuDUJRwwKVr6tQ6TkM3LnhmpCgwqe9qFh090JHJHEo4ML9RnzauFl23NqJFvhbNUUTTw$ 

%% Extract the data points to be used in the centroid calculation
% mag = []; usex = [];
% ind = 1;
% for ii = 1:length(peakRange)
%     if peakRange(ii) && ~any(usex==ii)
%         mag(ind) = y(ii);
%         usex(ind) = ii;
%         ind = ind+1;
%     end
%     if ii>1 && peakRange(ii) && ~any(usex==ii-1)
%         mag(ind) = y(ii-1);
%         usex(ind) = ii-1;
%         ind = ind+1;
%     end
%     if ii<length(y) && peakRange(ii) && ~any(usex==ii+1)
%         mag(ind) = y(ii+1);
%         usex(ind) = ii+1;
%         ind = ind+1;
%     end
% end
% [usex,si] = sort(usex);
% mag = mag(si);

mag = y(peakRange);
usex = find(peakRange == 1);

%% Work out the magnitude and angle of the sum of the vectors
% see https://urldefense.com/v3/__http://www.1728.org/vectutor.htm__;!!DZ3fjg!_lZLZrdVHoTKgDfZZNHO6tzdlLEuDUJRwwKVr6tQ6TkM3LnhmpCgwqe9qFh090JHJHEo4ML9RnzauFl23NqJFvhbNUUTTw$ 
xc = zeros(size(usex)); yc = xc;
for jj = 1:length(usex)
    xc(jj) = mag(jj)*cos(xRad(usex(jj)));
    yc(jj) = mag(jj)*sin(xRad(usex(jj)));
end

sumX = sum(xc);
sumY = sum(yc);

finalMag = sqrt(sumX^2+sumY^2);
finalAngle = atan2(sumY,sumX);
finalAngle = rad2deg(finalAngle);

% if finalAngle>90
%     finalAngle=-(finalAngle-90);
% elseif finalAngle<90
%     finalAngle=90-finalAngle;
% end

centroid = round(finalAngle,1);