% evaluate rate of strain tensor S(u)

% in 2D Cartesian coordinates:
% S(u) =  1/2* (grad(u)+grad(u)^T) 
%      =  1/2* ( 2*du/dx        du/dy + dv/dx )
%              ( du/dy + dv/dx  2*dv/dy       )
%      = (S11 S12)
%        (S21 S22)     


S11 = (1/2)* 2*(Su_ux*uh + ySu_ux);
S12 = (1/2)* (Su_uy*uh + ySu_uy + Sv_uy*vh + ySv_uy);
S21 = (1/2)* (Su_vx*uh + ySu_vx + Sv_vx*vh + ySv_vx);
S22 = (1/2)* 2*(Sv_vy*vh + ySv_vy);

% Note: S11 and S22 at xp,yp locations (pressure locations)
% S12, S21 at vorticity locations (x,y)

% for periodic boundary conditions S11(Npx+1,:)=S11(1,:)
% so S11 has size (Npx+1)*Npy; the last row are 'ghost' points equal to the
% first points. positions ([xp; xp(1)], yp)
% similarly, S22(:,Npy+1)=S22(:,1). positions (xp, [yp;yp(1)])

% 'cut-off' the double points in case of periodic BC
S11_p = reshape(S11,Nux_in+1,Nuy_in);
S11_p = S11_p(2:Nux_in+1,:);
S22_p = reshape(S22,Nux_in,Nuy_in+1);
S22_p = S22_p(:,2:Nuy_in+1);


% S12 is defined on the corners: size Nux_in*(Nuy_in+1), positions (xin,y)
% similarly S21 is size (Nux_in+1)*Nuy_in, positions (x,yin)
% get S12 and S21 at all corner points
S12_temp              = zeros(Nx+1,Ny+1);
S12_temp(1:Nx,:)      = reshape(S12,Nx,Ny+1);
S12_temp(Nx+1,:)      = S12_temp(1,:);

S21_temp              = zeros(Nx+1,Ny+1);
S21_temp(:,1:Ny)      = reshape(S21,Nx+1,Ny);
S21_temp(:,Ny+1)      = S21_temp(:,1);

% now interpolate to pressure points
S12_p = interp2(y',x,S12_temp,yp',xp); %interp2(x,y',S12_temp,xp,yp');
S21_p = interp2(y',x,S21_temp,yp',xp); %interp2(x,y',S21_temp,xp,yp');

% S21 and S12 should be equal!

% contour(xp,yp,S11');
% contour(xp,yp,S22');
% contour(xp,yp,S12');
% contour(xp,yp,S21');


%% invariants

% q = 1/2 * trace(S^2)
% r = -1/3 * trace(S^3) (= -det(S) only in 3D); in 2D we should get r=0

q = (1/2) * ( S11_p(:).^2 + S12_p(:).^2 + S21_p(:).^2 + S22_p(:).^2 );

% should be zero:
% r = (S11(:).^2+S12(:).*S21(:)).*S11(:) + (S11(:).*S12(:)+S12(:).*S22(:)).*S21(:) + ...
%     (S11(:).*S21(:)+S21(:).*S22(:)).*S12(:) + (S12(:).*S21(:)+S22(:).^2).*S22(:); %-(S11(:).*S22(:) - S12(:).*S21(:)); 

% figure
% contour(xp,yp,reshape(q,Npx,Npy)',25);
% figure
% contour(xp,yp,reshape(r,Npx,Npy)',25);


