%
fig_name = 'singular_values';

addpath('/Users/sanderse/Dropbox/work/Programming/libs/altmany-export_fig-cafc7c5/')

h = figure(1) % or whatever fig

%% general properties
set(gcf,'Color','w')
set(gca,'LineWidth',1)
set(gca,'FontSize',14)


%% lines
% get the axes handle
a = get(h, 'CurrentAxes')
% get the handles for children of the axes -- these are the data series handles
c = get(a, 'Children')

% change particular line:
% set(c(1),'Color','-');
set(c(3),'LineStyle','-');
set(c(3),'Marker','o');

set(c(1),'LineStyle','-');
set(c(1),'Marker','square');

% get color from colormap
cmap = get(gca,'ColorOrder');
set(c(1),'Color',cmap(2,:)); % pick 2nd color from map

%% restore color order
ax=gca;
ax.ColorOrderIndex=1; % or a different number depending where to start

%%
export_fig(fig_name,'-pdf','-transparent');