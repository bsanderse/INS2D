% post-processing
u    = reshape(uh,Nux_in,Nuy_in);
v    = reshape(vh,Nvx_in,Nvy_in);
pres = reshape(p,Npx,Npy);

line = {'r-','b-','m-','k-','g-','y-'};
color = char(line(j));

% we assume that Npx and Npy are even

% fix pressure jump (low Ct approximation)
deltap = 0.5*Ct;
% pres_ = pres - 0.5*(max(max(pres))+min(min(pres)));
% pres_ = (pres - 0.5*(max(max(pres))+min(min(pres))))/deltap;
pres_ = pres;

xr = xin(xmid);
yr = (yp(yrange(length(yrange)/2))+yp(yrange(length(yrange)/2+1)))/2;


figure(1)
% exact solution at the pressure points
% p_ex = - deltap / (2*pi) * ( atan((D/2-yr)./(xp-xr)) + atan((D/2+yr)./(xp-xr)) );
p_ex = -1 / (2*pi) * ( atan((D/2-yr)./(xp-xr)) + atan((D/2+yr)./(xp-xr)) );

% plot(xp,p_ex,'r-');
% hold on

px = (pres_(:,Npy/2+1)+pres_(:,floor(Npy/2)))/2;
plot(xp,px,color)
hold on
% legend('Exact','Present');
grid

%%
figure(2)
% exact solution at the pressure points
u_ex  = - deltap*p_ex - deltap*(xp>xr);
% plot(xp,u_ex,'r-');
% hold on
% ux = ((u(:,Nuy_in/2)+u(:,Nuy_in/2+1))/2 - u_inf)/deltap;
ux = (u(:,Nuy_in/2)+u(:,Nuy_in/2+1))/2;

plot(xin,ux,color)
hold on
uw = sqrt(1-Ct) - 1;
% duw = sqrt(1-Ct)-1;
% plot([min(x) max(x)],[uw uw]/deltap,'k--');



grid

%%
figure(3)
xd     = find(xin>=1,1); % point at which profile is plotted

disp(['wake velocity profile at: x=' num2str(xin(xd))]);
% u velocities are at grid lines
p_ex_y = -1/(2*pi) * (atan((D/2-yp)/xin(xd)) + atan((D/2+yp)/xin(xd)));
% u_ex_y =   - deltap*p_ex_y - deltap*(abs(yp)<D/2);
u_ex_y =   - p_ex_y - (abs(yp)<D/2);

% plot(yp,u_ex_y,'r-');
% hold on

% uy = (u(xd,:)-u_inf)/deltap;
uy = u(xd,:);
plot(yp,uy,color)
hold on

figure(4)
xd     = find(xp>=1,1); % point at which profile is plotted
% p is at pressure points
disp(['pressure profile at: x=' num2str(0.5*(xp(xd-1)+xp(xd)))]);

p_ex_y = -1/(2*pi) * (atan((D/2-yp)/xp(xd)) + atan((D/2+yp)/xp(xd)));
% plot(yp,p_ex_y,'r-');
% hold on
py = 0.5*(pres_(xd,:)+pres_(xd-1,:));

plot(yp,py,color);
hold on



figure(5)
vy = 0.5*(v(xd,:)+v(xd-1,:));
plot(yin,vy,color)
hold on

%%
% pp_AD_mesh_refinement;