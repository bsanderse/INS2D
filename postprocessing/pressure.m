
pres = reshape(p,Npx,Npy);

%BFS:
% set the pressure to zero at the corner of the step
% pres = pres-interp2(xp,yp,pres',0,0,'spline');

% shift = abs((pres(1,Npy/2)+pres(1,Npy/2+1))/2);
% shift = abs(pres(1,Npy/2));
% pres = pres + shift;

%LDC:
l = [0.3 0.17 0.12 0.11 0.09 0.07 0.05 0.02 0.0 -0.002];
if (floor(Nx/2)==Nx/2 && floor(Ny/2)==Ny/2)
    pres = pres-(pres(Nx/2+1,Ny/2+1)+pres(Nx/2,Ny/2))/2;
else
    pres = pres-pres(ceil(Nx/2),ceil(Ny/2));
end

% BFS:
% l = [0.01:0.01:0.1 0.12:0.02:0.24];

figure
% l = 25;
% l=min(p):1:max(p);
% l = linspace(-0.3,0.9,20);
contour(xp,yp,pres',l,'LineWidth',2);
axis equal
axis([x1 x2 y1 y2]);
xlabeltex('x',14);
ylabeltex('y',14);
grid
title('pressure');
colorbar
set(gca,'LineWidth',2)
