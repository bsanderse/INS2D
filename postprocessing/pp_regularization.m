% compute C(u,omega), the convection of omega with the velocity field u

% assume omega is known at x,y: omega_temp
vorticity;

%% interpolate omega_temp to u and v positions:
% omega_u  = interp2(x,y',omega_temp,xp,yp');

% to ux and uy positions
% omega_ux = zeros(Nx+1,Ny);
% omega_ux(1:Nx,1:Ny) = omega_u;
% omega_ux(Nx+1,:)    = omega_u(1,:);
% omega_uy = omega_temp(1:end-1,:);
% 
% % for v we can do the same:
% omega_vx = omega_temp(:,1:end-1);
% omega_vy = zeros(Nx,Ny+1);
% omega_vy(1:Nx,1:Ny) = omega_u;
% omega_vy(:,Ny+1)    = omega_u(:,1);

% to u and v positions
omega_u  = interp2(y',x,omega_temp,yp',xin); %interp2(x,y',omega_temp,xin,yp');
omega_v  = interp2(y',x,omega_temp,yin',xp); %interp2(x,y',omega_temp,xp,yin');


%%

uh_f = filter_convection(uh,Diffu_f,yDiffu_f,alfa); %uh + (alfa^2)*Re*(Diffu*uh + yDiffu);
vh_f = filter_convection(vh,Diffv_f,yDiffv_f,alfa);
cu   = uh_f;
cv   = vh_f;

omega_u_f = filter_convection(omega_u(:),Diffu_f,yDiffu_f,alfa);
omega_v_f = filter_convection(omega_v(:),Diffv_f,yDiffv_f,alfa);

% boundary conditions have to be changed if not periodic!!

w_ux   = Au_ux*omega_u_f+yAu_ux;                 % w at ux
uf_ux  = Iu_ux*cu+yIu_ux;                 % ubar at ux
duwdx  = Cux*(uf_ux.*w_ux);    

w_uy   = Au_uy*omega_u_f+yAu_uy;                 % w at uy
vf_uy  = Iv_uy*cv+yIv_uy;                 % vbar at uy
duwdy  = Cuy*(vf_uy.*w_uy);

w_vx   = Av_vx*omega_v_f+yAv_vx;                 % w at vx
uf_vx  = Iu_vx*cu+yIu_vx;                 % ubar at vx
dvwdx  = Cvx*(uf_vx.*w_vx);

w_vy   = Av_vy*omega_v_f+yAv_vy;                 % w at vy
vf_vy  = Iv_vy*cv+yIv_vy;                 % vbar at vy    
dvwdy  = Cvy*(vf_vy.*w_vy);


Cuw = filter_convection(duwdx + duwdy,Diffu_f,yDiffu_f,alfa);
Cvw = filter_convection(dvwdx + dvwdy,Diffv_f,yDiffv_f,alfa);

figure
labels=-0.15:0.015:0.15;
contour(xin,yp,reshape(Cuw,Nux_in,Nuy_in)',25)
caxis([-0.15 0.15]);
% contour(xp,yin,reshape(Cvw,Nvx_in,Nvy_in)',25)