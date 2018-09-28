% velocity values at pressure points

% restrict to velocities at pressure positions

up    = reshape( Bup*(Au_ux * uh + yAu_ux), Npx, Npy);
vp    = reshape( Bvp*(Av_vy * vh + yAv_vy), Npx, Npy);

qp    = sqrt(up.^2+vp.^2); 


skip = 1;
figure(1)
contour(xp,yp,qp',20);
hold on
quiver(xp(1:skip:end),yp(1:skip:end),up(1:skip:end,1:skip:end)',vp(1:skip:end,1:skip:end)',1);

axis equal
axis([x1 x2 y1 y2]);

xlabel('$x$','Interpreter','Latex');
ylabel('$y$','Interpreter','Latex');

colorbar

% airfoil
% hold on
% plot(x_cl,y_cl,'kx-')
% plot(x_k,y_k,'kx-');
% 
% 
hold off;
