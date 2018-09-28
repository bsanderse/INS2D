%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initial velocity and pressure

%    constant:    
        u  = ones(Nux_in,Nuy_in);
        v  = zeros(Nvx_in,Nvy_in);

    % BFS, extend inflow 
%         u  = (yu>=0).*(24*yu.*(1/2-yu));  
%         v  = zeros(Nvx_in,Nvy_in);

% same value as inflow:    
%         u     = kron(3*(yp>0)+1*(yp<=0),ones(Nux_in,1));  
%         v     = zeros(Nvx_in,Nvy_in);
%         v     = 0.01*sin(xv);

    % Eca:
    %     u   = erf(sigma*yu./xu);
    %     v   = 1/(sigma*sqrt(pi)) * (1 - exp(-(sigma*yv./xv).^2));

    % poiseuille like
%     u = kron(6*yp.*(1-yp),ones(Nux_in,1)) + ...
%         (16/(2*pi))*kron(2*yp-6*yp.^2+4*yp.^3,sin(2*pi*xin));
%     v = -16*kron(yin.^2-2*yin.^3+yin.^4,cos(2*pi*xp));

%       u = 6*yu.*(1-yu);
%       v = sin(pi*xv).^2.*cos(pi*2*yv);
%       u = 2*sin(2*pi*yu)*pi.*(16*xu.^4-(192/5)*xu.^5+32*xu.^6-(64/7)*xu.^7);
%       v = (64*xv.^3-192*xv.^4+192*xv.^5-64*xv.^6).*cos(2*pi*yv);
%         u = 0.2*sin(2*pi*yu).*sin(pi*xu).*cos(pi*xu) + 6*yu-6*yu.^2;
%         v = 0.2*(sin(pi*xv).^2-1/2).*cos(2*pi*yv);
%         u = 0.2*(1/2*cos(2*pi*xu).*sin(pi*yu)) + 6*yu-6*yu.^2;
%         v = 0.2*sin(2*pi*xv).*cos(pi*yv);
    % vortex, Knikker
%         C = 0.6944;
%         ru2 = xu.^2 + yu.^2;
%         rv2 = xv.^2 + yv.^2;
%         u = 1-yu*C.*exp(-0.5*ru2);%.*(ru2<L_x/4);
%         v = 1+xv*C.*exp(-0.5*rv2);%.*(rv2<L_y/4);
% %     %     

    % vortex a la Gresho
    % streamfunction psi=sin(pi*(x-x0)/L)^2 * sin(pi*(y-y0)/H)^2
    % take x0=y0=0, L=H=1 for a streamfunction centered at (0,0)
%     x0=0; y0=0; L=1; H=1;
%     u = -(2*pi/H)*(cos(pi*(xu-x0)/L).^2).*cos(pi*(yu-y0)/H).*sin(pi*(yu-y0)/H).*(xu<x0+L/2 & xu>x0-L/2 & yu<y0+H/2 & yu>y0-H/2);
%     v = +(2*pi/L)*cos(pi*(xv-x0)/L).*(cos(pi*(yv-y0)/H).^2).*sin(pi*(xv-x0)/L).*(xv<x0+L/2 & xv>x0-L/2 & yv<y0+H/2 & yv>y0-H/2);

    % sine pulse u
%           u  = (xu>L_x/4 & xu<3*L_x/4).*sin(pi*(xu-L_x/4)/(L_x/2)).^4;
%           u   = kron(ones(Nuy_in,1),u1D);
%           v   = zeros(Nvx_in,Nvy_in);

    % sine pulse v
    %       v   = (yv>L_y/4 & yv<3*L_y/4).*sin(pi*(yv-L_y/4)/(L_y/2)).^4;
    %       u   = kron(ones(Nuy_in,1),u1D);
    %       u   = zeros(Nux_in,Nuy_in);

%     double jet flow (domain [0,2*pi]x[0,2*pi])
%          d   = pi/15;
%          e   = 0.05;
%          u   = tanh( (yu-pi/2)/d) .* (yu<=pi) + tanh( (3*pi/2 - yu)/d) .* (yu>pi);
%          v   = e*sin(xv);
%          p   = zeros(Npx,Npy);
% 

    % decaying vortex (Taylor-Green)
%          u   = - sin(pi*xu).*cos(pi*yu);
%          v   = cos(pi*xv).*sin(pi*yv);
%          p   = (1/4)*(cos(2*pi*xpp) + cos(2*pi*ypp));

    % regularized LDC
%     u = 8*(xu.^4-2*xu.^3+xu.^2).*(4*yu.^3-2*yu);
%     v = -8*(4*xv.^3 - 6*xv.^2 + 2*xv).*(yv.^4-yv.^2);

    % adapted van Kan
%       u = (1/12)*(27/4)*xu.^2.*(3*xu.^2-8*xu+6)*6.*yu.*(1-yu);
%       v = -(27/4)*xv.*(xv-1).^2.*yv.^2.*(3-2*yv);
%         u = (128/45)*xu.^2.*(4*xu.^3-15*xu.^2+20*xu-10).*yu.*(yu-1);
%         v = -(256/27)*xv.*(xv.^3-3*xv.^2+3*xv-1).*yv.^2.*(2*yv-3);
%         xm = 1-2^(-1/3);
%         u = -(1/5)*xu.^2.*(15-40*xu+45*xu.^2-24*xu.^3+5*xu.^4).*yu.*(yu-1)/(xm*(xm-1)^4);
%         v = xv.*(xv-1).^4.*yv.^2.*(2*yv-3)/(xm*(xm-1)^4);

    % Kovasznay
%     lambda = Re/2-sqrt(Re^2/4+4*pi^2);
%     u = 1-exp(lambda*xu).*cos(2*pi*yu);
%     v = lambda/(2*pi) * exp(lambda*xv).*sin(2*pi*yv);

    % random on [-1,1]
%     u = -1 + 2*rand(Nux_in,Nuy_in);
%     v = -1 + 2*rand(Nvx_in,Nvy_in);

% for reversibility studies:
%      uh_start = u(:);
%      vh_start = v(:);

% dipole wall interaction (see e.g. Knikker)    
%     omega_e = 300;
%     r0 = 0.1;
%     xc1=0; yc1=0.1;
%     xc2=0; yc2=-0.1;
%     
%     u = 0.5*omega_e*( -(yu-yc1).*exp(-( (xu-xc1).^2 + (yu-yc1).^2 )/r0^2) + ...
%                       +(yu-yc2).*exp(-( (xu-xc2).^2 + (yu-yc2).^2 )/r0^2) );
%     v = 0.5*omega_e*(  (xv-xc1).*exp(-( (xv-xc1).^2 + (yv-yc1).^2 )/r0^2) + ...
%                       -(xv-xc2).*exp(-( (xv-xc2).^2 + (yv-yc2).^2 )/r0^2) );   
%     kin = 0.5*sum(Omu.*u(:).^2) + 0.5*sum(Omv.*v(:).^2)
% 
%     % scale velocity such that kinetic energy=2
%     u = u*sqrt(2)/sqrt(kin);
%     v = v*sqrt(2)/sqrt(kin);
%     omega_e = omega_e*sqrt(2)/sqrt(kin);
%     % this should equal 2:
%     kin     = 0.5*sum(Omu.*u(:).^2) + 0.5*sum(Omv.*v(:).^2);
% 
%     p = zeros(Np,1);


%%% initial turbulent quantities  

    % constant, same value as (constant) inlet:
%         kt  = kLe(1)*ones(Np,1);
%         e   = eLe(1)*ones(Np,1);

    % same value as inflow:    
    %    kt = kron(kLe,ones(Npx,1));                   
    %    e  = kron(eLe,ones(Npx,1));                
        
    % Eca:
    %     kt = kmax*(sigma_nu*ypp./xpp).^2 .* exp(1-(sigma_nu*ypp./xpp).^2);
    %     e  = 0.36*kmax^2/numax*exp(-(sigma_nu*ypp./xpp).^2);    

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%