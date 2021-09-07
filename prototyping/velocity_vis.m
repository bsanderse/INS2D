function velocity_vis(V,options)
xp = options.grid.xp;
yp = options.grid.yp;
[XX,YY] = meshgrid(xp,yp);


% [up,vp,qp] = get_velocity(V,t,options);
[up,vp] = get_inner_velocity(V,options);

% figure
% quiver(XX,YY,up',vp')

[xle,yle] = size(XX);
sparsity = 5;
x_s = 1:sparsity:xle;
y_s = 1:sparsity:yle;
XX_s = XX(x_s,y_s);
YY_s = YY(x_s,y_s);
U_s = up(y_s,x_s)';
V_s = vp(y_s,x_s)';
quiver(XX_s,YY_s,U_s,V_s)
hold on

% hori_x = linspace(xp(1),xp(end),length(xp)/sparsity);
% vert_y = linspace(yp(1),yp(end),length(yp)/sparsity);
hori_x = linspace(xp(1),xp(end),length(yp)/sparsity);
vert_y = linspace(yp(1),yp(end),length(yp)/sparsity);


x_ri = xp(1)*ones(size(vert_y));
x_le = xp(end)*ones(size(vert_y));
y_do = yp(1)*ones(size(hori_x));
y_up = yp(end)*ones(size(hori_x));

startx = [hori_x, x_ri, hori_x, x_le];
starty = [y_up, vert_y, y_do, vert_y];
streamline(XX,YY,up',vp',startx,starty);

% actuator
hold on
plot([2,2],[-0.5,0.5],'k-','LineWidth',3) % actuator disk
hold off

extra = xp(1)*sparsity;
x1 = xp(1)-extra;
x2 = xp(end)+extra;
y1 = yp(1)-extra;
y2 = yp(end)+extra;

axis equal
axis([x1 x2 y1 y2]);

% ylim([-2.025,2.025])