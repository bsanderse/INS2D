%% post-processing Couette results

line  = {'r-','b-','k-','m-','g-'};
color = char(line(j));


Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;
Npx    = options.grid.Npx;
Npy    = options.grid.Npy;    
Nu   = options.grid.Nu;
Nv   = options.grid.Nv;

uh   = V(1:Nu);
vh   = V(Nu+1:Nu+Nv);
u    = reshape(uh,Nux_in,Nuy_in);
v    = reshape(vh,Nvx_in,Nvy_in);
pres = reshape(p,Npx,Npy);

% shift pressure to get zero pressure in the centre
if (floor(Nx/2)==Nx/2 && floor(Ny/2)==Ny/2)
    pres_ = pres-(pres(Nx/2+1,Ny/2+1)+pres(Nx/2,Ny/2))/2;
else
    pres_ = pres-pres(ceil(Nx/2),ceil(Ny/2));
end

% streamfunction
psi = get_streamfunction(V,t,options);

% error in velocity field
% Couette solution:
u_exact = 1-2*yp';
u_error = max(abs(u(1,:) - u_exact))


%% figures
%
figure
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

figure
contour(x(2:end-1),y(2:end-1),reshape(psi,Nx-1,Ny-1)','k'); %labels,'LineWidth',1);


figure
l = 20; 
contour(xp,yp,pres_',l,'LineWidth',2);
axis equal
axis([x1 x2 y1 y2]);
xlabeltex('x',14);
ylabeltex('y',14);
grid
title('pressure');
colorbar
set(gca,'LineWidth',1)

figure
plot(time,umom-umom(1),'s-');
hold on
plot(time,vmom-vmom(1),'s-');
plot(time,k-k(1),'s-');
legend('u momentum','v momentum','kinetic energy');