% domain integration, assuming pressure right and lower/upper, inflow left


% momentum balance in y-direction
% steady flow, 
% pressure upper and lower balance (p=0)
u = reshape(uh,Nux_in,Nuy_in);
v = reshape(vh,Nvx_in,Nvy_in);
pres = reshape(p,Npx,Npy);


Bmap  = kron(speye(Nuy_in),BMx);
up    = reshape( Bmap*(Au_ux * uh + yAu_ux), Npx, Npy);
Bmap  = kron(BMy,speye(Nvx_in));
vp    = reshape( Bmap*(Av_vy * vh + yAv_vy), Npx, Npy);


Mom_x = - sum( -(uLe.^2).*hy) - sum((u(end,:).^2).*hy' ) + ...            % v^2 upper and lower contribution
        - sum( -up(:,1).*vp(:,1).*hx + up(:,end).*vp(:,end).*hx) + ...
        - sum( -pres(1,:).*hy' + pres(end,:).*hy')

% integrate over slightly smaller total volume
hx_new = hx;
hx_new(1) = 0.5*hx(1);
hx_new(end) = 0.5*hx(end);
Mom_y = - sum( -(v(:,1).^2).*hx_new + (v(:,end).^2).*hx_new ) + ...            % v^2 upper and lower contribution
        - sum( -up(1,:).*vp(1,:).*hy' + up(end,:).*vp(end,:).*hy') + ...   % u*v left and right side (at xp(1) and xp(end))
        - sum( -pres(:,1).*hx_new + pres(:,end).*hx_new)
