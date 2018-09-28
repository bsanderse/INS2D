% Kovasznay flow

% additional solve for pressure - unnecessary
% tn =0;
% A     = M*spdiags(1./Om,0,Nu+Nv,Nu+Nv)*G;
% [L,U] = lu(A);
% pressure_additional_solve;


p_ex = -0.5*exp(2*lambda*xpp);
Gpx_ex = -lambda*exp(2*lambda*xu);
Gpy_ex = zeros(Nv,1);
Gp_ex  = [Gpx_ex(:);Gpy_ex];

lambda = Re/2-sqrt(Re^2/4+4*pi^2);
u_ex = 1-exp(lambda*xu).*cos(2*pi*yu);
v_ex = lambda/(2*pi) * exp(lambda*xv).*sin(2*pi*yv);
V_ex = [u_ex(:);v_ex(:)];

diffu_ex = (1/Re)*(4*pi^2 - lambda^2) * exp(lambda*xu) .* cos(2*pi*yu);
diffv_ex = (1/Re)*(0.5*lambda^3/pi - 2*lambda*pi)* exp(lambda*xv).*sin(2*pi*yv);
diffuxx_ex = (1/Re)*(- lambda^2) * exp(lambda*xu) .* cos(2*pi*yu);

V = [uh;vh];

% figure
% contour(xp,yp,p_ex');

figure
surf(xp,yp,reshape(p_ex(:)-p,Npx,Npy)');
xlabel('x')
ylabel('y')

% figure
% surf(xp,yp,reshape(u_ex(:)-uh,Nux_in,Nuy_in)');
% xlabel('x')
% ylabel('y')

ep = p-p_ex(:);
errorpi(j) = max(abs(ep))
errorp2(j) = sqrt( sum(ep.^2)/Np)
errorpl(j) = sum( abs(ep) )/Np

eGp = Gp_ex-Om_inv.*(G*p+[y_px;y_py]);
errorGpi(j) = max(abs(eGp))
errorGp2(j) = sqrt( sum( eGp.^2)/(Nu+Nv) )

eV  = V-V_ex(:);
errorVi(j) = max(abs(eV))
errorV2(j) = sqrt( sum(eV.^2)/(Nu+Nv))
errorV1(j) = sum(abs(eV))/(Nu+Nv)

% uRi_i  = -lambda*exp(lambda*x2)*cos(2*pi*yp);
% ybc    = kron(uLe_i,Su_ux_BC.ybc1) + kron(uRi_i,Su_ux_BC.ybc2);
% ySu_ux = Su_ux_BC.Bbc*ybc;
dxx = Dux*( (1/Re)* Su_ux) * uh + Dux*( (1/Re)* ySu_ux);

ediffu = diffu_ex(:) - Omu_inv.*(Diffu*uh + yDiffu);
ediffv = diffv_ex(:) - Omv_inv.*(Diffv*vh + yDiffv);
ediffuxx = diffuxx_ex(:) - Omu_inv.*dxx;

% leave out last point
temp = reshape(ediffuxx,Nux_in,Nuy_in);
ediffuxx = temp(1:Nux_in-1,:);

errorDui(j) = max(abs(ediffu))
errorDvi(j) = max(abs(ediffv))
errorDu2(j) = sqrt( sum(ediffu.^2)/Nu )
errorDv2(j) = sqrt( sum(ediffv.^2)/Nv )
errorDuxxi(j) = max(abs(ediffuxx(:)))
errorDuxx2(j) = sqrt( sum(ediffuxx(:).^2)/((Nux_in-1)*Nuy_in) )