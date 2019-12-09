%% real-time plotting lid-driven cavity 

line  = {'r-','b-','k-','m-','g-'};
color = char(line(j));


Npx    = options.grid.Npx;
Npy    = options.grid.Npy;    

figure(1)
[up,vp,qp] = get_velocity(V,t,options);
qp = reshape(qp,Npx,Npy);
labels= linspace(-1,1,20);
% labels=20;
contour(xp,yp,qp',labels);
hold on
quiver(xp,yp,up',vp',1);
axis square
colorbar
grid
hold off