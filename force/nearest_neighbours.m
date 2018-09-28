function [ xn, yn, neighbour, weight ] = nearest_neighbours( xp, yp, x, y )
%NEAREST_NEIGHBOURS find grid points xn,yn nearest to a point xp,yp 
% given a Cartesian mesh x, y

% number of points
n = 4;

xn = zeros(2,1);
yn = zeros(2,1);

% x-direction
indx = find(x-xp>=0,1);

% if (abs(x(indx)-xp)<eps)
%     xn(1:2) = xp; % co-inciding with grid line
%     n  = n - 1;
%     dx1 = 0;
%     dx2 = 1;
% else    
    xn  = x(indx-1:indx);
    dx2 = abs(x(indx)-xp);
    dx1 = abs(x(indx-1)-xp);
% end

% y-direction
indy = find(y-yp>=0,1);

if (abs(y(indx)-yp)<eps)
    yn(1:2) = yp; % co-inciding with grid line
    n  = n - 1;
    dy1 = 0;
    dy2 = 1;    
else    
    yn = y(indy-1:indy);
    dy2 = abs(y(indy)-yp);
    dy1 = abs(y(indy-1)-yp);
end

% xn = [xn
% neighbour=0;
neighbour.x(1:4) = [xn(1); xn(2); xn(1); xn(2)];
neighbour.y(1:4) = [yn(1); yn(1); yn(2); yn(2)];
weight(1) = dx2*dy2;
weight(2) = dx1*dy2;
weight(3) = dx1*dy1;
weight(4) = dx2*dy1;
weight    = weight/((dx1+dx2)*(dy1+dy2))

end

