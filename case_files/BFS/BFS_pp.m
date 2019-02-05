% post-processing backward facing step

Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;
Npx    = options.grid.Npx;
Npy    = options.grid.Npy;
yu     = options.grid.yu;

uh   = V(1:Nu);
vh   = V(Nu+1:Nu+Nv);
u    = reshape(uh,Nux_in,Nuy_in);
v    = reshape(vh,Nvx_in,Nvy_in);


%% calculate length of reattachment zones
figure(1)

start = 8; % prevent first points to be taken
plot(x(2:end),u(:,1),'b');
hold on
plot(x(2:end),u(:,end),'r');
grid
legend(['u at y=' num2str(yu(1,1))],['u at y=' num2str(yu(1,end))])
xlabel('x'); ylabel('u')

ubottom = u(start:end,1);

% find first point that has positive velocity
il = find(min(ubottom,0)==0,1);
if (~isempty(il))
    % interpolate to obtain intersection point
    X1(j) = x(il+start-1) + (-ubottom(il-1)/(ubottom(il)-ubottom(il-1)))*(x(il+start)-x(il+start-1));
    disp(['X1: ' num2str(X1*2)]);
    
    plot(X1,0,'ko','MarkerSize',10)
end

utop = u(:,end);
il = find(max(utop,0)==0,1);
if (~isempty(il))
    X2(j) = x(il) + (-utop(il-1)/(utop(il)-utop(il-1)))*(x(il+1)-x(il));
    disp(['X2: ' num2str(X2*2)]);
    plot(X2,0,'ko','MarkerSize',10)
end
il = find(min(utop(il:end),0)==0,1) + il - 1;
if (~isempty(il))
    X3(j) = x(il) + (-utop(il-1)/(utop(il)-utop(il-1)))*(x(il+1)-x(il));
    disp(['X3: ' num2str(X3*2)]);
    plot(X3,0,'ko','MarkerSize',10)
end

%% profiles at x=7 and x=15
figure
hold on

i7 = find(x==7);
if (isempty(i7))
    disp('x=7 is not a point in the results');
else
    plot(u(i7-1,:),yp,'kx-');
end

i15=find(x==15);
if (isempty(i15))
    disp('x=15 is not a point in the results');
else   
    plot(u(i15-1,:),yp,'bx-');
end

y_Gartling = [0.5 0.45 0.4 0.35 0.3 0.25 0.20 0.15 0.10 0.05 0.0 -0.05 -0.10 -0.15 -0.2 -0.25 -0.3 -0.35 -0.4 -0.45 -0.5];
u_Gartling7 = [0.0 -0.038 -0.049 -0.032 0.015 0.092 0.204 0.349 0.522 0.709 0.885 1.024 1.105 1.118 1.062 0.948 0.792 0.613 0.428 0.232 0.000];
u_Gartling15 = [0.0 0.101 0.202 0.304 0.408 0.512 0.613 0.704 0.779 0.831 0.853 0.844 0.804 0.737 0.649 0.547 0.438 0.328 0.218 0.109 0.00];

% u_ex15 = (a/8)*(yp+L_y/2).*(L_y/2-yp);

plot(u_Gartling7,y_Gartling,'ro-');
plot(u_Gartling15,y_Gartling,'go-');
legend('x=7, current','x=15, current','x=7, Gartling', 'x=15, Gartling');


%% velocity
[up,vp,qp] = get_velocity(V,t,options);

% BFS:
% l = [0.05 0.1 0.15 0.2 0.4 0.6 0.8 1 1.2 1.4];
list = 0.1:0.05:1.3;
% list = 25;

figure
contour(xp,yp,qp',list)
axis([x1 x2 y1 y2]);
axis equal
colorbar
set(gca,'LineWidth',2)
grid
title('velocity');
colorbar

%% pressure
pres = reshape(p,Npx,Npy);

%BFS:
% set the pressure to zero at the corner of the step
pres  = pres-interp2(xp,yp,pres',0,0,'spline');
shift = abs((pres(1,Npy/2)+pres(1,Npy/2+1))/2);
% shift = abs(pres(1,Npy/2));
pres = pres + shift;

% BFS:
l = [0.01:0.01:0.1 0.12:0.02:0.24];

figure
contour(xp,yp,pres',l,'LineWidth',2);
axis equal
axis([x1 x2 y1 y2]);
xlabeltex('x',14);
ylabeltex('y',14);
grid
title('pressure');
colorbar
set(gca,'LineWidth',2)

%% vorticity
omega = get_vorticity(V,t,options);
omega = reshape(omega,Nx-1,Ny-1);

figure
labels = -8:2:10;
contour(x(2:end-1),y(2:end-1),reshape(omega,Nx-1,Ny-1)',labels,'LineWidth',1);
% 
axis equal
axis([x1 x2 y1 y2]);
% 
xlabeltex('x',14);
ylabeltex('y',14);
colorbar
grid
title('vorticity')
set(gca,'LineWidth',2);
% hold off;

%% streamfunction
psi = get_streamfunction(V,t,options);
labels = 25;
contour(x(2:end-1),y(2:end-1),reshape(psi,Nx-1,Ny-1)',labels,'LineWidth',2);
axis equal
axis([x1 x2 y1 y2]);
xlabel('x');
ylabel('y');
colorbar
grid
title('streamfunction')
set(gca,'LineWidth',2)

