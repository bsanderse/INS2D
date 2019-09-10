%% post-processing lightly loaded actuator disk simulations
% comparison with analytical results
% see also section 5.3 in
% ECNS: Energy-Conserving Navier-Stokes Solver - Verification of steady laminar flows
% B. Sanderse, ECN-E–11-042, June 2011

% note: we assume that Npx and Npy are even
if (rem(Npx,2)~=0 || rem(Npy,2)~=0)
    error('number of volumes in x- and y-direction should be even for postprocessing');
end

%% get grid variables
Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;
Nu     = options.grid.Nu;
Nv     = options.grid.Nv;

xin  = options.grid.xin;
yin  = options.grid.yin;

%% get solution variables
uh   = V(1:Nu);
vh   = V(Nu+1:end);
u    = reshape(uh,Nux_in,Nuy_in);
v    = reshape(vh,Nvx_in,Nvy_in);
pres = reshape(p,Npx,Npy);
[up,vp,qp] = get_velocity(V,t,options);

figure
contour(xp,yp,qp');

figure
contour(xp,yp,pres');


%% plot settings
line = {'r-','b-','m-','k-','g-','y-'};
color = char(line(j));

%% actuator settings
x_c  = options.force.x_c;
y_c  = options.force.y_c;

% pressure jump (low Ct approximation)
deltap = 0.5*Ct;

% pres_ = pres - 0.5*(max(max(pres))+min(min(pres)));
% pres_ = (pres - 0.5*(max(max(pres))+min(min(pres))))/deltap;
pres_ = pres;

%% actuator position on grid
xmid = find(xin>=x_c,1);
xr   = xin(xmid);

if (rem(Nuy_in,2)==0)
    % Nuy_in even
    ymid   = Nuy_in/2;
    % number of volumes in disk (assuming locally uniform in y-dir)
    nd     = round(D/hy(ymid));
    yrange = ymid-floor(nd/2)+1:ymid+floor(nd/2);
else
    % Nuy_in odd
    ymid   = floor(Nuy_in/2)+1;
    % number of volumes in disk (assuming locally uniform in y-dir)
    nd     = D/hy(ymid);
    yrange = ymid-floor(nd/2):ymid+floor(nd/2);
end

yr = (yp(yrange(length(yrange)/2))+yp(yrange(length(yrange)/2+1)))/2;


%% compare analytical and numerical solutions

%% centerline pressure

% analytical solution
p_ex = -deltap / (2*pi) * ( atan((D/2-yr)./(xp-xr)) + atan((D/2+yr)./(xp-xr)) );
% numerical solution
px   = (pres_(:,Npy/2+1)+pres_(:,floor(Npy/2)))/2;


figure(1)
plot(xp,p_ex);
hold on

plot(xp,px,color)
hold on
% legend('Exact','Present');
grid
title('pressure at centerline');

%% centerline velocity

% analytical solution
u_ex = 1 - p_ex - deltap*(xp>xr);
% numerical solution
ux   = (u(:,Nuy_in/2)+u(:,Nuy_in/2+1))/2;

figure(2)
plot(xp,u_ex);
hold on
plot(xin,ux,color)
% hold on
uw = sqrt(1-Ct) - 1;
% duw = sqrt(1-Ct)-1;
% plot([min(x) max(x)],[uw uw]/deltap,'k--');

title('velocity at centerline');
grid


%% wake profile velocity

figure(3)

xwake  = 1;
xd     = find(xin>=xwake,1); % point at which profile is plotted

disp(['wake velocity profile at: x=' num2str(xin(xd))]);

% u velocities are at grid lines
p_ex_y = -deltap/(2*pi) * (atan((D/2-yp)/xin(xd)) + atan((D/2+yp)/xin(xd)));
% u_ex_y =   - deltap*p_ex_y - deltap*(abs(yp)<D/2);
u_ex_y = 1 - p_ex_y - deltap*(abs(yp)<D/2);

plot(yp,u_ex_y);
hold on

% uy = (u(xd,:)-u_inf)/deltap;
uy = u(xd,:);
plot(yp,uy,color)
hold on
title('velocity in wake');
grid


%% wake profile pressure

figure(4)

xwake  = 1;
xd     = find(xp>=xwake,1); % point at which profile is plotted
% p is at pressure points
disp(['pressure profile at: x=' num2str(0.5*(xp(xd-1)+xp(xd)))]);

p_ex_y = -deltap/(2*pi) * (atan((D/2-yp)/xp(xd)) + atan((D/2+yp)/xp(xd)));
plot(yp,p_ex_y);
hold on
py = 0.5*(pres_(xd,:)+pres_(xd-1,:));

plot(yp,py,color);
hold on

title('pressure in wake');
grid

%% vertical velocity profile 

figure(5)
vy = 0.5*(v(xd,:)+v(xd-1,:));
plot(yin,vy,color)
hold on

%%
% pp_AD_mesh_refinement;