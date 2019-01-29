function psi = get_streamfunction(V,t,options)
% compute streamfunction from a Poisson equation nabla^2 psi = -omega

global uBC vBC

BC = options.BC;
Nu     = options.grid.Nu;
Nux_in = options.grid.Nux_in;
Nvx_in = options.grid.Nvx_in;
Nx     = options.grid.Nx;
Ny     = options.grid.Ny;

hx = options.grid.hx;
hy = options.grid.hy;
x  = options.grid.x;
y  = options.grid.y;
xp = options.grid.xp;
yp = options.grid.yp;

Wv_vx = options.discretization.Wv_vx;
Wu_uy = options.discretization.Wu_uy;

uh = V(1:Nu);
vh = V(Nu+1:end);

% boundary values by integrating around domain
% start with psi=0 at lower left corner

% u = d psi / dy; integrate low->up
if (strcmp(BC.u.left,'dir')) 
%     u1 = interp1(y,uLe,yp);
    u1 = uBC(x(1),yp,t,options);
elseif (strcmp(BC.u.left,'pres') || strcmp(BC.u.left,'per'))
    u1 = uh(1:Nux_in:end);
end
psiLe  = cumsum(hy.*u1);
psiUpLe= psiLe(end);
psiLe  = psiLe(1:end-1);

% v = -d psi / dx; integrate left->right
if (strcmp(BC.v.up,'dir')) 
%     v1 = interp1(x,vUp,xp);
    v1 = vBC(xp,y(end),t,options);
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
    u2 = uBC(x(end),yp,t,options);
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
    v2 = vBC(xp,y(1),t,options);
elseif (strcmp(BC.v.low,'pres') || strcmp(BC.v.low,'per'))
    v2 = vh(1:Nvx_in);
end
psiLo  = psiLoRi + cumsum(hx(end:-1:1).*v2(end:-1:1)); 
psiLoLe= psiLo(end); 
psiLo  = psiLo(end-1:-1:1);

if (abs(psiLoLe)>1e-12)
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

omega = get_vorticity(V,t,options);

% solve streamfunction from Poisson equaton
psi  = -Apsi\(omega + yApsi);