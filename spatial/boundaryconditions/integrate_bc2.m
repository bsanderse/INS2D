% integrate BC
% used for RK methods

% this only works for boundary conditions which are exponential in time

% Taylor vortex
% fac1 = (-2*pi^2/Re);
% if (i_RK==1)
%     fac2 = 1 + dt*c_RK(1)*fac1;
% elseif (i_RK==2)
%     fac2 = 1 + dt*c_RK(2)*fac1 + (dt^2)*A_RK(1,1)*A_RK(2,2)*fac1^2;
% elseif (i_RK==3)
%     fac2 = 1 + dt*c_RK(3)*fac1 + ...
%               (dt^2)*(A_RK(3,2)*A_RK(1,1) + A_RK(3,3)*(A_RK(2,1)+A_RK(2,2)))*fac1^2 + ...
%               (dt^3)*(A_RK(3,3)*A_RK(2,2)*A_RK(1,1))*fac1^3;
% end

% van Kan
if (tn>0)
    
    der1 = 1/(tn^2);
    der2 = 1/(tn^4) - 2/(tn^3);
    der3 = 1/(tn^6) - 6/(tn^5) + 6/(tn^4);

    if (i_RK==1)
        fac2 = 1 + dt*c_RK(1)*der1;
    elseif (i_RK==2)
        fac2 = 1 + dt*c_RK(2)*der1 + (dt^2)*A_RK(1,1)*A_RK(2,2)*der2;
    elseif (i_RK==3)
        fac2 = 1 + dt*c_RK(3)*der1 + ...
                  (dt^2)*(A_RK(3,2)*A_RK(1,1) + A_RK(3,3)*(A_RK(2,1)+A_RK(2,2)))*der2 + ...
                  (dt^3)*(A_RK(3,3)*A_RK(2,2)*A_RK(1,1))*der3;
    end

else 
    
    fac2 = 1;
    
end

uLo = uLo_n*fac2;
uUp = uUp_n*fac2;
uLe = uLe_n*fac2;
uRi = uRi_n*fac2;

uLo_i = uLo_i_n*fac2;
uUp_i = uUp_i_n*fac2;
uLe_i = uLe_i_n*fac2;
uRi_i = uRi_i_n*fac2;

vLo = vLo_n*fac2;
vUp = vUp_n*fac2;
vLe = vLe_n*fac2;
vRi = vRi_n*fac2;

vLo_i = vLo_i_n*fac2;
vUp_i = vUp_i_n*fac2;
vLe_i = vLe_i_n*fac2;
vRi_i = vRi_i_n*fac2;
    
   