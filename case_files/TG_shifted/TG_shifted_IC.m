function [u,v,p,options] = TG_shifted_IC(t,options)
% initial velocity field Taylor-Green

xu = options.grid.xu;
yu = options.grid.yu;
xv = options.grid.xv;
yv = options.grid.yv;
xpp = options.grid.xpp;
ypp = options.grid.ypp;

u   = - sin(pi*xu).*cos(pi*yu);
v   = cos(pi*xv).*sin(pi*yv);

% pressure: should in principle NOT be prescribed. will be calculated if
% p_initial=1
p   = (1/4)*(cos(2*pi*xpp) + cos(2*pi*ypp));
    
end