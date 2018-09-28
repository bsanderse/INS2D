% post-processing 2D results
u    = reshape(uh,Nux_in,Nuy_in);
v    = reshape(vh,Nvx_in,Nvy_in);
pres = reshape(p,Npx,Npy);

line = {'r-','b-','m-','k-','g-','y-'};
color = char(line(j));

% we assume that Npx and Npy are even

% fix pressure jump (low Ct approximation)
deltap = 0.5*Ct;

%% case 1
% construct a line normal to AD
alfa   = 30;
xwake1 = linspace(x1,x2,46)';
ywake1 = tand(alfa)*xwake1 -1 + tand(alfa);
xint   = x1+cosd(alfa)^2;
yint   = y1+cosd(alfa)*sind(alfa);
swake1 = sign(xwake1-xint).*sqrt((xwake1-xint).^2 + (ywake1-yint).^2);

px = interp2(xp,yp,pres',xwake1,ywake1);
px2 = pres(:,5);
figure
plot(swake1,px,'rx-')
hold on
plot(xp,px2,'mx-')

Vx = interp2(xp,yp,qp',xwake1,ywake1);
figure
plot(swake1,Vx,'rx-')

%% case 2
figure
plot(xp,pres(:,1))