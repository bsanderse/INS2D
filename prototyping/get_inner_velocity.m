function [up,vp] = get_inner_velocity(V,options)
% hard-coded for actuator_unsteady_ROM

xp = length(options.grid.xp);
yp = length(options.grid.yp);

up_1 = V(1:xp*yp);
vp_1 = V(xp*yp+1:end);

% up_2 = reshape(up_1,yp,xp);
% vp_2 = reshape(vp_1,yp+1,xp);
up_2 = reshape(up_1,xp,yp)';
vp_2 = reshape(vp_1,xp,yp+1)';

up = [up_2(:,1) (up_2(:,1:end-1)+up_2(:,2:end))/2];
vp = (vp_2(1:end-1,:)+vp_2(2:end,:))/2;

up = up';
vp = vp';


