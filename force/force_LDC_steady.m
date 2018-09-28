% point force

Shih_f     = xv.^4 - 2*xv.^3 + xv.^2;
Shih_f1    = 4*xv.^3 - 6*xv.^2 + 2*xv;
Shih_f2    = 12*xv.^2 - 12*xv + 2;
Shih_f3    = 24*xv - 12;
Shih_g     = yv.^4 - yv.^2;
Shih_g1    = 4*yv.^3 - 2*yv;
Shih_g2    = 12*yv.^2 - 2;

Shih_F     = 0.2*xv.^5 - 0.5*xv.^4 + (1/3)*xv.^3;
Shih_F1    = -4*xv.^6 + 12*xv.^5 - 14*xv.^4 + 8*xv.^3 - 2*xv.^2;
Shih_F2    = 0.5*(Shih_f.^2);
Shih_G1    = -24*yv.^5 + 8*yv.^3 - 4*yv;

Fy         = (-8/Re)*(24*Shih_F + 2*Shih_f1.*Shih_g2 + Shih_f3.*Shih_g) - ...
             64*(Shih_F2.*Shih_G1 - Shih_g.*Shih_g1.*Shih_F1);

% this works for both 2nd and 4th order method
Fy         = -Omv.*Fy(:);

Fx         = zeros(Nu,1);

         
% at pressure points, for pressure solution         
Shih_p_f   = xpp.^4 - 2*xpp.^3 + xpp.^2;
Shih_p_f1  = 4*xpp.^3 - 6*xpp.^2 + 2*xpp;
Shih_p_f2  = 12*xpp.^2 - 12*xpp + 2;
Shih_p_f3  = 24*xpp - 12;
Shih_p_g   = ypp.^4 - ypp.^2;
Shih_p_g1  = 4*ypp.^3 - 2*ypp;
Shih_p_g2  = 12*ypp.^2 - 2;
Shih_p_g3  = 24*ypp;


Shih_p_F   = 0.2*xpp.^5 - 0.5*xpp.^4 + (1/3)*xpp.^3;
Shih_p_F1  = -4*xpp.^6 + 12*xpp.^5 - 14*xpp.^4 + 8*xpp.^3 - 2*xpp.^2;
Shih_p_F2  = 0.5*(Shih_p_f.^2);
Shih_p_G1  = -24*ypp.^5 + 8*ypp.^3 - 4*ypp;   
         
%%
% reasonable plot for Re=100:
% Fy = reshape(Omv_inv.*Fy,Nvx_in,Nvy_in);
% contour(xp,yin,Fy',-2.5:0.1:0.5)

