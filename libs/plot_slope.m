function plot_slope(x1,y1,x2,slope,triangle,orientation)
% plot_slope(x1,y1,x2,slope,triangle,orientation)
% plot slope in log-log plot
% triangle=1 draws triangle instead of line only
% orientation=0 means lower side horizontal, orientation=1 means upper side
% horizontal

if (nargin<4)
    error('wrong inputs');
end
if (nargin==4)
    triangle=0;
end
if (nargin==4 || nargin==5)
    orientation=0;
end

if (triangle~=0 && triangle~=1)
    warning('triangle specification invalid');
end

if (orientation~=0 && orientation~=1)
    warning('orientation invalid');
end
y2 = exp(slope*(log(x2)-log(x1)) + log(y1));

hold on

loglog([x1 x2],[y1 y2],'k-');

if (triangle==1)
    if (orientation==0)
        loglog([x1 x2],[y1 y1],'k-');
        loglog([x2 x2],[y1 y2],'k-');
        text(exp((log(x1)+log(x2))/2),exp((2*log(y1)+log(y2))/3),num2str(slope));
    end
    if (orientation==1)
        loglog([x1 x1],[y1 y2],'k-');
        loglog([x1 x2],[y2 y2],'k-');
        text(exp((3*log(x1)+log(x2))/4),exp((log(y1)+2*log(y2))/3),num2str(slope));

    end
end