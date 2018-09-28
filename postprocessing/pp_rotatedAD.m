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
pres_ = (pres - 0.5*(max(max(pres))+min(min(pres))))/deltap;
% pres_ = pres;

% xr = xin(xmid);
% yr = (yp(yrange(length(yrange)/2))+yp(yrange(length(yrange)/2+1)))/2;

% coordinates through centerline wake
xwake1 = linspace(x1,x2,100)';
ywake1 = linspace(-0.5*(x2-x1)*tand(alfa),0.5*(x2-x1)*tand(alfa),100)';
swake1 = sign(ywake1).*sqrt(xwake1.^2 + ywake1.^2);
% ywake1 = ywake1 + y_c;

% coordinates 1D behind turbine
xm = D*cosd(alfa);
ym = D*sind(alfa); 
xm1 = xm-(y2-ym)*tand(alfa); xm2 = min(x2,xm+(ym-y1)*tand(alfa));
ym1 = y2; ym2 = max(y1,ym-(x2-xm)/tand(alfa));

xwake2 = linspace(xm1,xm2,100)';
ywake2 = linspace(ym1,ym2,100)';
swake2 = sign(ywake2-ym).*sqrt((xwake2-xm).^2 + (ywake2-ym).^2);

%%
figure(1)
% exact solution at the pressure points
% p_ex = - deltap / (2*pi) * ( atan((D/2-yr)./(xp-xr)) + atan((D/2+yr)./(xp-xr)) );
p_ex = -1 / (2*pi) * ( atan((D/2-y_c)./(swake1-x_c)) + atan((D/2+y_c)./(swake1-x_c)) );

plot(swake1,p_ex,'k-');
hold on

% px = (pres_(:,Npy/2+1)+pres_(:,floor(Npy/2)))/2;
px = interp2(xp,yp,pres_',xwake1,ywake1);
plot(swake1,px,color)
legend('Exact','Present');
grid

%%
figure(2)
% exact solution at the pressure points
u_ex  = 1- deltap*p_ex - deltap*(swake1>0);
ux = interp2(xin,yp,u',xwake1,ywake1);
vx = interp2(xp,yin,v',xwake1,ywake1);
% Vx = interp2(xp,yp,qp',xwake1,ywake1);
Vt = ux*cosd(alfa) + vx*sind(alfa);

plot(swake1,Vt,'rx-')
hold on
plot(swake1,u_ex,'kx-')

uw = sqrt(1-Ct) - 1;
% duw = sqrt(1-Ct)-1;
% plot([min(x) max(x)],[uw uw]/deltap,'k--');
grid

%%
figure(3)
p_ex_y = -1/(2*pi) * (atan((D/2-swake2)/D) + atan((D/2+swake2)/D));
% u_ex_y =   - deltap*p_ex_y - deltap*(abs(yp)<D/2);
u_ex_y =  1 - deltap*p_ex_y - deltap*(abs(swake2)<D/2);

uy = interp2(xin,yp,u',xwake2,ywake2);
vy = interp2(xp,yin,v',xwake2,ywake2);
Vt = uy*cosd(alfa) + vy*sind(alfa);

plot(swake2,u_ex_y,'kx-');
hold on
plot(swake2,Vt,'rx-')

%%
figure(4)
py = interp2(xp,yp,pres_',xwake2,ywake2);

plot(swake2,p_ex_y,'kx-');
hold on
plot(swake2,py,'rx-')



% figure(5)
% Vn = -uy*sind(alfa) + vy*cosd(alfa);
% plot(swake2,Vn,color)
% hold on

%%
% pp_AD_mesh_refinement;