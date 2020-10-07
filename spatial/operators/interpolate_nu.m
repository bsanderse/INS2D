function [nu_t_ux,nu_t_uy,nu_t_vx,nu_t_vy] = interpolate_nu(nu_t,options)
% interpolate the scalar field nu_t at pressure locations (xp,yp)
% to locations needed in computing the diffusive terms, i.e. the u_x, u_y,
% v_x and v_y locations

%% option 1: using the operators defined in operator_interpolate_viscosity

Anu_ux   = options.discretization.Anu_ux;
yAnu_ux  = options.discretization.yAnu_ux;
nu_t_ux  = Anu_ux * nu_t + yAnu_ux;

Anu_uy   = options.discretization.Anu_uy;
yAnu_uy  = options.discretization.yAnu_uy;
nu_t_uy  = Anu_uy * nu_t + yAnu_uy;

Anu_vx   = options.discretization.Anu_vx;
yAnu_vx  = options.discretization.yAnu_vx;
nu_t_vx  = Anu_vx * nu_t + yAnu_vx;

Anu_vy   = options.discretization.Anu_vy;
yAnu_vy  = options.discretization.yAnu_vy;
nu_t_vy  = Anu_vy * nu_t + yAnu_vy;

%% option 2: construct interpolation operators here 

% Npx = options.grid.Npx;
% Npy = options.grid.Npy;        
% Nux_in = options.grid.Nux_in;
% Nuy_in = options.grid.Nuy_in;
% Nvx_in = options.grid.Nvx_in;
% Nvy_in = options.grid.Nvy_in;        
% 
% xp = options.grid.xp;
% yp = options.grid.yp;
% xin = options.grid.xin;
% yin = options.grid.yin;
% x  = options.grid.x;
% y  = options.grid.y;
% hx = options.grid.hx;
% hy = options.grid.hy;
% 
% BC = options.BC;

% %% ux; corresponding to S11
% if (strcmp(BC.u.left,'per') && strcmp(BC.u.right,'per'))
%     % ux positions where we need nu:
%     % x_out = ([xp(1)-0.5*(hx(1)+hx(end)); xp], yp)     
%     % y_out = yp;
%     % so nu at xp,yp is almost correct, except for the additional first
%     % column
%     % so we add extra points for periodic BC
%     nu_t_ux               = zeros(Nux_in+1,Nuy_in);
%     nu_t_ux(2:Nux_in+1,:) = nu_t;
%     nu_t_ux(1,:)          = nu_t(end,:); 
%     
% elseif (strcmp(BC.u.left,'dir') && strcmp(BC.u.right,'pres'))
%     % ux positions are correct
%     % add extra points for outflow BC, which are such that when taking
%     % derivatives, a zero gradient results
%     nu_t_ux               = zeros(Nux_in+1,Nuy_in);
%     nu_t_ux(1:Nux_in,:)   = nu_t;
%     nu_t_ux(Nux_in+1,:)   = nu_t(end,:);
%         
% else
%     error('BC not implemented in strain_tensor.m');
% end
% 
% %% uy; corresponding to S12
% 
% % construct a larger matrix of (Npx+2)*(Npy+2) which can be used for interpolation
% nu_t_ghost            = zeros(Npx+2,Npy+2);
% nu_t_ghost(2:Npx+1,2:Npy+1) = nu_t; % set interior
% nu_t_ghost(2:Npx+1,1)   = nu_t(:,end);
% nu_t_ghost(2:Npx+1,end) = nu_t(:,1);
% nu_t_ghost(1,2:Npy+1)   = nu_t(end,:);
% nu_t_ghost(end,2:Npy+1) = nu_t(1,:);
% nu_t_ghost(1,1)         = nu_t(end,end);
% nu_t_ghost(end,end)     = nu_t(1,1);
% nu_t_ghost(1,end)       = nu_t(end,1);
% nu_t_ghost(end,1)       = nu_t(1,end);
% 
% if (strcmp(BC.u.left,'per') && strcmp(BC.u.right,'per'))
%     
%     xp_in                  = [xp(1)-hx(1);xp;xp(end)+hx(end)];
%     yp_in                  = [yp(1)-hy(1);yp;yp(end)+hy(end)];    
%     nu_t_uy                = interp2(yp_in',xp_in,nu_t_ghost,y',xin);
%     
% elseif (strcmp(BC.u.left,'dir') && strcmp(BC.u.right,'pres'))
% 
% else
%     error('BC not implemented in strain_tensor.m');
% end
% 
% 
% %% vx; corresponding to S21
% 
% if (strcmp(BC.v.low,'per') && strcmp(BC.v.up,'per'))
%     
%     xp_in                  = [xp(1)-hx(1);xp;xp(end)+hx(end)];
%     yp_in                  = [yp(1)-hy(1);yp;yp(end)+hy(end)];       
%     nu_t_vx                = interp2(yp_in',xp_in,nu_t_ghost,yin',x);
% 
% elseif (strcmp(BC.v.low,'pres') && strcmp(BC.v.up,'pres'))
% 
%     
% else
%     error('BC not implemented in strain_tensor.m');
% end
% 
% 
% %% vy; corresponding to S22
% 
% if (strcmp(BC.v.low,'per') && strcmp(BC.v.up,'per'))
%     % vy positions where we need nu:
%     % x_out = xp;    
%     % y_out = (xp,[yp(1)-0.5*(hy(1)+hy(end)); yp])     
%     % so nu at xp,yp is almost correct, except for the additional first
%     % column
%     % so we add extra points for periodic BC    
%     nu_t_vy               = zeros(Nvx_in,Nvy_in+1);
%     nu_t_vy(:,2:Nvy_in+1) = nu_t;
%     nu_t_vy(:,1)          = nu_t(:,end);
%     
% elseif (strcmp(BC.v.low,'pres') && strcmp(BC.v.up,'pres'))
%     % vy positions are correct
%     % add extra points for outflow BC
%     % these are copied from the interior, which leads to a zero gradient,
%     % which is fine as a no-stress boundary condition is used here    
%     nu_t_vy               = zeros(Nvx_in,Nvy_in+1);
%     nu_t_vy(:,2:Nvy_in)   = nu_t;
%     nu_t_vy(:,1)          = nu_t(:,1);
%     nu_t_vy(:,Nvy_in+1)   = nu_t(:,end); 
%     
%     
% else
%     error('BC not implemented in strain_tensor.m');
% end






