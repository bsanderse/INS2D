% figure(111)
% hold on
% x = linspace(0,100);
% y = x.^2;
% 
% plot(x,y)
% 
% axes('position',[.65,.65,.35,.25])
% box on
% ioi1 = (x<70) & (x>20);
% ioi2 = (y<40) & (y>10);
% plot(ioi1,ioi2)
% 
% figure(666)
% ax1 = axes('Position',[.1 .1 .7 .7]);
% ax2 = axes('Position',[.2 .4 .25 .25]);
% 
% plot(ax1,x,y)
% plot(ax2,x(40:50),y(40:50))

close all
figure(45)
x = linspace(0,100);
y = x.^2;
                        ax1 = axes('Position',[.1 .1 .12 .8]);
                        ax2 = axes('Position',[.35 .1 .6 .8]);
% ax1 = axes('Position',[.1 .1 .12 .8]);
% ax2 = axes('Position',[.3 .1 .65 .8]);
plot(ax1,x,y)
plot(ax2,x,y)
xlim(ax1,[14,25.5])
ylabel(ax1,'hahh')
ylabel(ax2,'hahh')

% xlim(ax2,[0,4])