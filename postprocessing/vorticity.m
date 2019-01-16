% vorticity

omega = get_vorticity(V,t,options);
% 
% % velocity at inner points
% diagpos = 0;
% if (strcmp(BC.u.left,'pres') )
%     diagpos = 1;
% end
% if (strcmp(BC.u.right,'per') && strcmp(BC.u.left,'per') ) % like pressure left
%     diagpos = 1;
% end
% 
% B1D    = spdiags(ones(Nx-1,1),diagpos,Nx-1,Nux_in);
% B2D    = kron(speye(Ny),B1D);
% 
% uh_in  = B2D*uh;
% 
% 
% diagpos = 0;
% if (strcmp(BC.v.low,'pres') )
%     diagpos = 1;
% end
% if (strcmp(BC.v.low,'per') && strcmp(BC.v.up,'per') ) % like pressure low
%     diagpos = 1;
% end
% 
% B1D    = spdiags(ones(Ny-1,1),diagpos,Ny-1,Nvy_in);
% B2D    = kron(B1D,speye(Nx));
% 
% vh_in  = B2D*vh;
% 
% 
% omega  = Wv_vx*vh_in - Wu_uy*uh_in;


figure
% labels = -5:1:5;
% LDC:
labels = [-5 -4 -3 -2 -1 -0.5 0 0.5 1 2 3];
% labels = 30;

% BFS:
% labels = -8:2:10;

% labels = 25;
% labels=min(omega(:)):20:max(omega(:));
% labels = -omega_e:20:omega_e;
% labels(16)=[];
% labels = -0.01:0.001:0.01;
% labels=-40:1:40;
% labels(41)=[];
% labels = 20;
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
set(gca,'LineWidth',1);
% hold off;


% %% for periodic BC (e.g. doublejet flow)
% if (strcmp(BC.u.left,'per') && strcmp(BC.v.low,'per'))
% 
% % point value
% omega_p = Wv_vx_p*vh - Wu_uy_p*uh;
% % conserved quantity
% omega   = Wv_vx*vh - Wu_uy*uh;
% 
% if (order4==1)
%     % second order conserved quantity
%     omega1  = Wv_vx1*vh - Wu_uy1*uh;
% end
%     
% % labels=-4:0.5:4;
% labels=25;
% 
% omega_temp = zeros(Nx+1,Ny+1);
% omega_temp(1:Nx,1:Ny) = reshape(omega_p,Nx,Ny);
% omega_temp(Nx+1,1:Ny) = omega_temp(1,1:Ny);
% omega_temp(1:Nx,Ny+1) = omega_temp(1:Nx,1);
% omega_temp(Nx+1,Ny+1) = omega_temp(1,1);
% % contour(x(1:end-1),y(1:end-1),reshape(omega,Nx,Ny)',labels);
% contour(x,y,reshape(omega_temp,Nx+1,Ny+1)',labels);
% % 
% % axis equal
% % axis([x1 x2 y1 y2]);
% % fontsize=16;
% % xlabeltex('$x$',fontsize);
% % ylabeltex('$y$',fontsize);
% % set(gca,'fontsize',fontsize-3);
% % set(gca,'fontsize',fontsize-3);
% end