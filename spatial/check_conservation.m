function [maxdiv,umom,vmom,k] = check_conservation(V,t,options)
% check mass, momentum and energy conservation properties of velocity field
Nu = options.grid.Nu;
Nv = options.grid.Nv;

uh = V(1:Nu);
vh = V(Nu+1:Nu+Nv);

global uBC vBC;

BC  = options.BC;
if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(t,options);
end

M   = options.discretization.M;
yM  = options.discretization.yM;

Omu = options.grid.Omu;
Omv = options.grid.Omv;
x   = options.grid.x;
y   = options.grid.y;
xp  = options.grid.xp;
yp  = options.grid.yp;
hx  = options.grid.hx;
hy  = options.grid.hy;
gx  = options.grid.gx;
gy  = options.grid.gy;

uLe_i = uBC(x(1),yp,t,options);
uRi_i = uBC(x(end),yp,t,options);
vLo_i = vBC(xp,y(1),t,options);
vUp_i = vBC(xp,y(end),t,options);


%% check if new velocity field is divergence free (mass conservation)
maxdiv  = max(abs( M*V + yM ));


%% calculate total momentum
umom    = sum(Omu.*uh);
vmom    = sum(Omv.*vh);

% add boundary contributions in case of Dirichlet BC
if (strcmp(BC.u.left,'dir'))
    umom = umom + sum(uLe_i.*hy)*gx(1);
    % 4th order
%     umom(n) = umom(n) + sum(uLe_i.*(alfa*hy*gx(1)-hy3*(gx(1)+gx(2))));
end
if (strcmp(BC.u.right,'dir'))
    umom = umom + sum(uRi_i.*hy)*gx(end);
end
if (strcmp(BC.v.low,'dir'))
    vmom = vmom + sum(vLo_i.*hx)*gy(1);
end
if (strcmp(BC.v.up,'dir'))
    vmom = vmom + sum(vUp_i.*hx)*gy(end);
end


%% calculate total kinetic energy
k     = 0.5*sum(Omu.*uh.^2) + 0.5*sum(Omv.*vh.^2); % this equals 0.5*(V')*(Omega.*V);

% add boundary contributions in case of Dirichlet BC
if (strcmp(BC.u.left,'dir'))
    k = k + 0.5*sum((uLe_i.^2).*hy)*gx(1);
end
if (strcmp(BC.u.right,'dir'))
    k = k + 0.5*sum((uRi_i.^2).*hy)*gx(end);
end
if (strcmp(BC.v.low,'dir'))
    k = k + 0.5*sum((vLo_i.^2).*hx)*gy(1);
end
if (strcmp(BC.v.up,'dir'))
    k = k + 0.5*sum((vUp_i.^2).*hx)*gy(end);
end

% if (strcmp(BC.u.left,'per') && strcmp(BC.v.low,'per'))
% 
%     convection;
%     u_ux  = reshape(u_ux,Nx+1,Ny);
%     v_vy  = reshape(v_vy,Nx,Ny+1);
%     if (order4==0)
%         % point value of kinetic energy
%         k_p   = 0.5*u_ux(2:end,:).^2 + 0.5*v_vy(:,2:end).^2;
%         % volume integrated quantity
%         k2(n) = sum(sum(k_p(:).*Omp));
%     elseif (order4==1)
%         u_ux3 = reshape(u_ux3,Nx+3,Ny);
%         v_vy3 = reshape(v_vy3,Nx,Ny+3);
%         % point value of kinetic energy    
%         k_p   = 0.5*(beta*u_ux(2:end,:) + (1-beta)*u_ux3(3:end-1,:)).^2 + ...
%                 0.5*(beta*v_vy(:,2:end) + (1-beta)*v_vy3(:,3:end-1)).^2; 
%         % fourth order approximation to integral value
%         % reuse operator for vorticity
%         k2(n) = sum(Wint*(Omp.*k_p(:)));
%     end
% 
% 
%     % vorticity and enstrophy
%     vorticity;
% 
%     % the conserved quantity; in the second order case this is also the volume integrated quantity
%     omega_total(n)  = sum(omega);
% 
%     omega_total3(n) = sum(omega(1:end/2));
% 
%     if (order4==1)
%     % conserved quantity for second order method
%     omega_total1(n) = sum(omega1);
%     % volume quantity such that fourth order and total area=L_x*L_y
% 
%     omega_total2(n) = sum(Wint*(Omvort.*omega_p));
%     omega_total3(n) = sum(Wint(1:end/2,:)*(Omvort.*omega_p));
%     end
%     % enstrophy(n)   = 0.5*sum(Omvort.*(omega.^2)); 
% 
%     % for periodic BC:
%     enstrophy(n)  = -0.5*(uh'*(Re)*Diffu*uh + vh'*(Re)*Diffv*vh);
% 
% end

% 
% % reversibility error, inviscid flow
% % if (n<=nt/2+1)
%     uh_save(n,1:Nu) = uh;
%     vh_save(n,1:Nv) = vh;
% if (n>nt/2+1)
%     n_rev = nt+2-n;
% %     rev_error(n) = sqrt(sum(Omu.*(uh - uh_save(n_rev,:)').^2)/sum(Omu)) + sqrt(sum(Omv.*(vh - vh_save(n_rev,:)').^2)/sum(Omv));
%     rev_error(n) = sqrt(sum(Omu.*(uh - uh_save(n_rev,:)').^2) + sum(Omv.*(vh - vh_save(n_rev,:)').^2)/(2*k(n)));
% % elseif (n==nt+1)
% % %     rev_error(n) = sum(Omu.*(uh - uh_start).^2)/sum(Omu) + sum(Omv.*(vh - vh_start).^2)/sum(Omv);
% %     rev_error(n) = sqrt(sum(Omu.*(uh - uh_start).^2) + sum(Omv.*(vh - vh_start').^2)/(2*k(n)));
% end

end
