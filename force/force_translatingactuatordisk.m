% if (exist('n','var')==1 )

% actuator disk
Ct = 0.5;   % thrust coefficient
D  = 1;     % diameter
deltap = -0.5*Ct;  % pressure jump over AD
 
% period   = 2; % period


% previous position
if (i_RK == 1)
    
    if ( n==2 ) % first time step
        dx_old= 1; % moving to the right
        xc1_old = x_c; xc2_old = x_c;
        yc1_old = y_c-D/2; yc2_old = y_c+D/2;    
    else
        dx_old  = dx;
        xc1_old = xc1; yc1_old = yc1;
        xc2_old = xc2; yc2_old = yc2;
    end
    
end

% current position
xc1 = x_c + sin(t);
yc1 = y_c - D/2;
xc2 = x_c + sin(t);
yc2 = y_c + D/2;


dx = sign(cos(t)); % moving to right (+) or left (-)

% outer indices of boundary to limit computations in force_integral
imin = find(x>min([xc1_old;xc1]),1)-2;
imax = find(x>max([xc1_old;xc1]),1)+1;
jmin = find(y>yc1,1)-2;
jmax = find(y>yc2,1)+1;

% when movement reverses the boundaries should be wider:
if (dx~=dx_old)
    if (dx>0)
       imin = find(x>x_c-1,1)-2;
    else
       imax = find(x>x_c+1,1)+1;
    end
end
       
bnd  = [imin,imax,jmin,jmax];

if (dx_old == dx) % same sign

    % is dx-dx_old is small enough, the following polygon suffices
    % note we subtract the center location to get the right definition for
    % r
    % CCW orientation of polygon
  
    % moving to the right
    AD   = [[xc1_old; xc1; xc2; xc2_old; xc1_old]-x_c,...
            [yc1_old; yc1; yc2; yc2_old; yc1_old]-y_c];
    if (dx<0)
        % moving to the left is obtained as a simple permutation of moving
        % to the right
        identity = speye(size(AD,1));
        perm = identity(end:-1:1,:); % permutation matrix
        AD   = [perm*AD(:,1) perm*AD(:,2)];
    end

    
    % compute intersectional area of disk with underlying mesh
    % midpoints of mesh for u-component is given by (xin,yp); the FV boundaries
    % are given by ([ ],y) in case of dirichlet/pressure BC
    Area_u = force_integral(AD,[xp;xp(end)+hx(end)]-x_c,y-y_c,bnd,0);
    
    % difference with exact area
    S_ex = abs(polygonArea(AD));
    if (abs(sum(sum(Area_u))-S_ex)>1e-10);
        warning('total impulse not conserved');        
    end
    
    % take into account normal velocity leads to
    Area_u = Area_u*(c_RK(i_RK)*dt)/abs(xc1-xc1_old);

else
    % angle that splits the area in two
    x_max = -dx + x_c; % -1 or 1
    
    if (dx_old>0) %right side
        % area 1; CCW orientation of polygon
        AD1 = [[xc1_old; x_max; x_max; xc2_old; xc1_old]-x_c,...
               [yc1; yc1; yc2; yc2; yc1]-y_c];
        % area 2; keep CCW orientation of polygon
        AD2 = [[xc1; x_max; x_max; xc2; xc1]-x_c,...
               [yc1; yc1; yc2; yc2; yc1]-y_c];
    else % left side
        % area 1; CCW orientation of polygon
        AD1 = [[x_max; xc1_old; xc2_old; x_max; x_max]-x_c,...
               [yc1; yc1; yc2; yc2; yc1]-y_c];
        % area 2; keep CCW orientation of polygon
        AD2 = [[x_max; xc1; xc2; x_max; x_max]-x_c,...
               [yc1; yc1; yc2; yc2; yc1]-y_c];           
    end
            
    Area_u1 = force_integral(AD1,[xp;xp(end)+hx(end)]-x_c,y-y_c,bnd,0);
    Area_u2 = force_integral(AD2,[xp;xp(end)+hx(end)]-x_c,y-y_c,bnd,0);
    
    % difference with exact area
    S_ex1 = abs(polygonArea(AD1));
    S_ex2 = abs(polygonArea(AD2));
    if (abs(sum(sum(Area_u1))-S_ex1)>1e-10 || abs(sum(sum(Area_u2))-S_ex2)>1e-10)
        warning('total impulse not conserved');     
        keyboard;
    end
    
    % take into account normal velocity leads to
    Area_u1 = Area_u1*(c_RK(i_RK)*dt)/abs(x_max-xc1_old);
    Area_u2 = Area_u2*(c_RK(i_RK)*dt)/abs(x_max-xc1);
    
    Area_u  = Area_u1 + Area_u2;
end


Fx_moving = deltap*Area_u;
Fx_moving = Fx_moving(:);
