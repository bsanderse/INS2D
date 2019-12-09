function layout(figno,fontsize,xtext,ytext,ztext)
% layout figure
h=figure(figno);
set(gca,'fontsize',fontsize);

% get the axes handle
% a = get(h, 'CurrentAxes')
% % get the handles for children of the axes -- these are the data series handles
% c = get(a, 'Children')
% 
% set(c,'fontsize',fontsize)

xlabeltex(xtext,fontsize+2)
ylabeltex(ytext,fontsize+2)
if (nargin==5)
    zlabeltex(ztext,fontsize+2)
end

