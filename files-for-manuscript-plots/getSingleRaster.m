a = findobj(gcf,'type','axes');
s = a(1);
a = findobj(s,'type','line');

x = a(end).XData;
y = a(end).YData;

% target 1 spikes
t1_inds = find(y(1:3:end)+0.5 <= 10);
t1_inds = 3*(t1_inds-1)+1;
t1_inds = sort([t1_inds,t1_inds+1,t1_inds+2]);
x1 = x(t1_inds); y1 = y(t1_inds);

% target 2 spikes
t2_inds = find(y(1:3:end)+0.5 > 10);
t2_inds = 3*(t2_inds-1)+1;
t2_inds = sort([t2_inds,t2_inds+1,t2_inds+2]);
x2 = x(t2_inds); y2 = y(t2_inds)-10;

figure; plot(x1/10000,y1,'k'); ylim([0.5 10.5]); xlim([0.450 0.750])
figure; plot(x2/10000,y2,'k'); ylim([0.5 10.5]); xlim([0.450 0.750])
