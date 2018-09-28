% streamfunction

% boundary values by integrating around domain
% start with psi=0 at lower left corner

% u = d psi / dy; integrate low->up
if (strcmp(BC.u.left,'dir')) 
%     u1 = interp1(y,uLe,yp);
    u1 = uBC(x(1),yp,t,Re);
elseif (strcmp(BC.u.left,'pres') || strcmp(BC.u.left,'per'))
    u1 = uh(1:Nux_in:end);
end
psiLe  = cumsum(hy.*u1);
psiUpLe= psiLe(end);
psiLe  = psiLe(1:end-1);

% v = -d psi / dx; integrate left->right
if (strcmp(BC.v.up,'dir')) 
%     v1 = interp1(x,vUp,xp);
    v1 = vBC(xp,y(end),t,Re);
elseif (strcmp(BC.v.up,'pres'))
    v1 = vh(end-Nvx_in+1:end);
elseif (strcmp(BC.v.up,'per'))
    v1 = vh(1:Nvx_in);
end
psiUp  = psiUpLe - cumsum(hx.*v1); 
psiUpRi= psiUp(end);
psiUp  = psiUp(1:end-1);

% u = d psi / dy; integrate up->lo
if (strcmp(BC.u.right,'dir')) 
%     u2 = interp1(y,uRi,yp);
    u2 = uBC(x(end),yp,t,Re);
elseif (strcmp(BC.u.right,'pres'))
    u2 = uh(Nux_in:Nux_in:end);
elseif (strcmp(BC.u.right,'per'))
    u2 = uh(1:Nux_in:end);
end
psiRi  = psiUpRi - cumsum(hy(end:-1:1).*u2(end:-1:1));
psiLoRi= psiRi(end);
psiRi  = psiRi(end-1:-1:1);

% v = -d psi / dx; integrate right->left
if (strcmp(BC.v.low,'dir')) 
%     v2 = interp1(x,vLo,xp);
    v2 = vBC(xp,y(1),t,Re);
elseif (strcmp(BC.v.low,'pres') || strcmp(BC.v.low,'per'))
    v2 = vh(1:Nvx_in);
end
psiLo  = psiLoRi + cumsum(hx(end:-1:1).*v2(end:-1:1)); 
psiLoLe= psiLo(end); 
psiLo  = psiLo(end-1:-1:1);

if (abs(psiLoLe)>1e-14)
    disp(['warning: contour integration of psi not consistent: ' num2str(abs(psiLoLe))]);
end

% solve del^2 psi = -omega
% only dirichlet boundary conditions because we calculate streamfunction at
% inner points only

% x-direction
diag1 = (1./hx).*ones(Nx,1);
Q1D   = spdiags([-diag1 diag1],[0 1],Nx,Nx+1);
Qx_BC = BC_general(Nx+1,Nx-1,2, ...
                                      'dir','dir',hx(1),hx(end));
% extend to 2D
Q2Dx  = kron(speye(Ny-1),Q1D*Qx_BC.B1D);
yQx   = kron(speye(Ny-1),Q1D*Qx_BC.Btemp)*...
         (kron(psiLe.*ones(Ny-1,1),Qx_BC.ybc1)+kron(psiRi.*ones(Ny-1,1),Qx_BC.ybc2));

     
% y-direction
diag1 = (1./hy).*ones(Ny,1);
Q1D   = spdiags([-diag1 diag1],[0 1],Ny,Ny+1);
Qy_BC = BC_general(Ny+1,Ny-1,2, ...
                                      'dir','dir',hy(1),hy(end));
% extend to 2D
Q2Dy  = kron(Q1D*Qy_BC.B1D,speye(Nx-1));
yQy   = kron(Q1D*Qy_BC.Btemp,speye(Nx-1))*...
         (kron(Qy_BC.ybc1,psiLo.*ones(Nx-1,1))+kron(Qy_BC.ybc2,psiUp.*ones(Nx-1,1)));
     
     
Apsi  = Wv_vx*Q2Dx + Wu_uy*Q2Dy;
yApsi = Wv_vx*yQx + Wu_uy*yQy;

if (exist('omega','var')==0)
    vorticity;
end
    

psi  = -Apsi\(omega + yApsi);


figure(4)
% Re=1000:
labels = [1.5e-3 1e-3 5e-4 2.5e-4 1e-4 5e-5 1e-5 1e-6 0 -1e-10 -1e-5 -1e-4 -1e-2 -3e-2 -5e-2 -7e-2 -9e-2 -0.1 -0.11 -0.115 -0.1175];
% Re=5000:
% labels = [-0.1175 -0.115 -0.11 -0.1 -0.09 -0.07 -0.05 -0.03 -0.01 -1e-4 -1e-5 -1e-7 0 1e-6 1e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 1.5e-3 3e-3];
% Re=10000:
% labels = [-0.1175 -0.115 -0.11 -0.1 -0.09 -0.07 -0.05 -0.03 -0.01 -1e-4 -1e-5 -1e-7 0 1e-6 1e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 1.5e-3 3e-3];

% BFS, Re=800
% labels = [-0.03 -0.025 -0.02 -0.015 -0.01 -0.005 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.49 0.5 0.502 0.504];
% labels = 60;
% labels =-0.12:0.005:0;
% labels1 = 1-(0:0.025:1).^2;
% labels2 = 1+(0:0.025:1).^2;
% labels = [flipud(labels1(2:end)');labels2'];
% labels=0.45:0.002:0.54;
contour(x(2:end-1),y(2:end-1),reshape(psi,Nx-1,Ny-1)',labels);

axis equal
axis([x1 x2 y1 y2]);

xlabel('x');
ylabel('y');

% actuator disks
% hold on
% contour(xu,yu,reshape(Fx,Nux_in,Nuy_in),min(Fx),'Color','black');
% h1=plot([x(xmid+1) x(xmid+1)],[yp(yrange(1)) yp(yrange(length(yrange)))],'k-');
% set(h1,'LineWidth',1.5);
% h2=plot([x(2*xmid+1) x(2*xmid+1)],[yp(yrange(1)+floor(nd/2)) yp(yrange(end)+floor(nd/2))],'k-');
% set(h2,'LineWidth',1.5);

% point sources
% hold on
% theta=0:0.01:2*pi;
% plot(x1+L_x/2+deltax*cos(theta),y1+L_y/2+deltax*sin(theta),'k-');
% 
% hold off
