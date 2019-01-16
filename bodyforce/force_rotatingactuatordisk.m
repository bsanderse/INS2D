% if (exist('n','var')==1 )

% actuator disk
Ct = 0.1;   % thrust coefficient
D  = 1;     % diameter
deltap = 0.5*Ct;  % pressure jump over AD (minus sign is taken into account in integrals)
alfa_max = (pi/3); % maximum angle
period   = 2; % period
freq     = 2*pi/period; % rotational frequency


% previous position
if (i_RK == 1)
    
    if ( n==2 ) % first time step
        alfa_old = 0;
        dalfa_old= 1;
        xc1_old  = x_c; xc2_old = x_c;
        yc1_old  = y_c+D/2; yc2_old = y_c-D/2;    
    else
        alfa_old = alfa; dalfa_old = dalfa;
        xc1_old = xc1; yc1_old = yc1;
        xc2_old = xc1; yc2_old = yc1;
    end
    
end

% current position
% angle (given by triangle wave in time)
period   = period/2;
alfa  = -alfa_max*(2/period)*(t-period*floor(t/period+1/2)).*(-1).^(floor(t/period-1/2));
% increasing or decreasing angle
dalfa = sign(cos(pi*t/period)); %alfa_max*(2/per)*dt;

if (abs(alfa-alfa_old)>2*alfa_max)
    error('time step larger than total turning angle');
end

% upper end
xc1 = x_c - (D/2)*sin(alfa);
yc1 = y_c + (D/2)*cos(alfa);
% lower end
xc2 = x_c + (D/2)*sin(alfa);
yc2 = y_c - (D/2)*cos(alfa);


% outer indices of boundary to limit computations in force_integral
imin = find(x>min([x_c;xc1_old;xc1]),1)-2;
imax = find(x>max([x_c;xc1_old;xc1]),1)+1;
jmin = find(y>min([y_c;yc1_old;yc1]),1)-2;
jmax = find(y>max([y_c;yc1_old;yc1]),1)+1;
bnd  = [imin,imax,jmin,jmax];

if (dalfa_old == dalfa) % same sign

    % is delta_alfa=alfa-alfa_old is small enough, the following polygon suffices
    % note we subtract the center location to get the right definition for
    % r
    % CCW orientation of polygon
    if (dalfa>0) %increasing angle
        AD   = [[x_c; xc1_old; xc1; x_c]-x_c,[y_c; yc1_old; yc1; y_c]-y_c];
    else
        AD   = [[x_c; xc1; xc1_old; x_c]-x_c,[y_c; yc1; yc1_old; y_c]-y_c];
    end

    
    % compute intersectional area of disk with underlying mesh
    % midpoints of mesh for u-component is given by (xin,yp); the FV boundaries
    % are given by ([ ],y) in case of dirichlet/pressure BC
%     Area_u   = comp_Area(xAD,yAD,[xp;xp(end)+hx(end)],y,1);

    % exact integral int int -sin(theta) dr dtheta
    B_ex = (D/2)*(cos(pi/2+alfa)-cos(pi/2+alfa_old));
    % exact integral int int cos(theta) dr dtheta
    C_ex = (D/2)*(sin(pi/2+alfa)-sin(pi/2+alfa_old));
   
    Area_u = force_integral(AD,[xp;xp(end)+hx(end)]-x_c,y-y_c,bnd,2);
    Area_v = force_integral(AD,x-x_c,[yp(1)-hy(1);yp;yp(end)+hy(end)]-y_c,bnd,3);
    
    % difference larger than 10%
    if (abs((sum(sum(Area_u))-B_ex)/B_ex)>0.1 || abs((sum(sum(Area_v))-C_ex)/C_ex)>0.1) 
        warning('total impulse not conserved');
    end

else
    % angle that splits the area in two
    alfa_mid = sign(alfa)*alfa_max;
    
    xc1_mid = x_c - (D/2)*sin(alfa_mid);
    yc1_mid = y_c + (D/2)*cos(alfa_mid);  
    
    if (dalfa>0) %increasing angle
        % area 1; CCW orientation of polygon
        AD1 = [[x_c; xc1_old; xc1_mid; x_c]-x_c,[y_c; yc1_old; yc1_mid; y_c]-y_c];
        % area 2; keep CCW orientation of polygon
        AD2 = [[x_c; xc1; xc1_mid; x_c]-x_c,[y_c; yc1; yc1_mid; y_c]-y_c];
    else
        % area 1; CCW orientation of polygon
        AD1 = [[x_c; xc1_mid; xc1_old; x_c]-x_c,[y_c; yc1_mid; yc1_old; y_c]-y_c];
        % area 2; keep CCW orientation of polygon
        AD2 = [[x_c; xc1_mid; xc1; x_c]-x_c,[y_c; yc1_mid; yc1; y_c]-y_c];
    end        
        
    Area_u1 = force_integral(AD1,[xp;xp(end)+hx(end)]-x_c,y-y_c,bnd,2);
    Area_v1 = force_integral(AD1,x-x_c,[yp(1)-hy(1);yp;yp(end)+hy(end)]-y_c,bnd,3);
    Area_u2 = force_integral(AD2,[xp;xp(end)+hx(end)]-x_c,y-y_c,bnd,2);
    Area_v2 = force_integral(AD2,x-x_c,[yp(1)-hy(1);yp;yp(end)+hy(end)]-y_c,bnd,3);

    Area_u  = Area_u1 + Area_u2;
    Area_v  = Area_v1 + Area_v2;
end

% mirror to get lower side; 
% assume that: (1) AD is in the middle of the domain
%              (2) Nx, Ny are even
%              (3) dirichlet/pressure for u, and pressure/pressure for v

Area_dum = zeros(Nux_in,Nuy_in);
Area_dum(1:end-1,1:end) = Area_u(end-1:-1:1,end:-1:1);
Area_u = Area_u + Area_dum;

Area_dum = zeros(Nvx_in,Nvy_in);
Area_dum(1:end,1:end) = Area_v(end:-1:1,end:-1:1);
Area_v = Area_v + Area_dum;

Fx_moving = (deltap/freq)*Area_u;
Fy_moving = (deltap/freq)*Area_v;

Fx_moving = Fx_moving(:);
Fy_moving = Fy_moving(:);

