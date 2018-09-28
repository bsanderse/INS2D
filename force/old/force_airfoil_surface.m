% prescribe Cp and Cf along airfoil surface

% currently naca only
airfoil_type = '2412';
aoa          = 0;   % in angle of attack
Re           = 1e7; % Reynolds number or 0 for inviscid
xfoil        = 1;

% load Cp and Cf from xfoil or own panel method
load_Cp_Cf;
    

% coordinates of leading edge
x_c       = L_x/2;
y_c       = L_y/2;

x_k       = x_k + x_c;
y_k       = y_k + y_c;

% find intersections of airfoil geometry with mesh.
% airfoil_intersections;
[xi, yi, panel, param] = geometry_intersection(x,y,x_k,y_k,Sk);
figure(1)
plot(x_k,y_k,'rx-')
plot(xi,yi,'bo-')

% for each finite volume intersected by a line segment find the integrated force
ni    = length(xi);
fx_as = zeros(Nx,Ny);
fy_as = zeros(Nx,Ny);

CpSx  = Cp.*Sk.*vecn(:,1);
CpSy  = Cp.*Sk.*vecn(:,2);
CfSx  = Cf.*Sk.*vect(:,1);
CfSy  = Cf.*Sk.*vect(:,2);

% a negative Cp is in the direction of the unit vector
Cx    = -CpSx+CfSx;
Cy    = -CpSy+CfSy;

for i=1:ni-1
    k1 = panel(i);
    k2 = panel(i+1);
    
    % find indices of finite volume
    % points inside the FV
    
        xtemp = 0.5*(xi(i)+xi(i+1));
        ytemp = 0.5*(yi(i)+yi(i+1)); 

        % closest pressure point=index FV
        [valx indx] = min(abs(xp-xtemp));
        [valy indy] = min(abs(yp-ytemp));

    % find integrated force
    if (k2-k1==0)
      fx_as(indx,indy) = Cx(k1)*(param(i+1)-param(i));
      fy_as(indx,indy) = Cy(k1)*(param(i+1)-param(i));
    elseif (k2-k1==1)
      fx_as(indx,indy) = Cx(k1)*(1-param(i)) + Cx(k2)*param(i+1);
      fy_as(indx,indy) = Cy(k1)*(1-param(i)) + Cy(k2)*param(i+1);
    elseif (k2-k1>1)
      fx_as(indx,indy) = Cx(k1)*(1-param(i)) + Cx(k2)*param(i+1) + ...
                         sum(Cx(k1+1:k2-1));
      fy_as(indx,indy) = Cy(k1)*(1-param(i)) + Cy(k2)*param(i+1) + ...
                         sum(Cy(k1+1:k2-1));            
    else
        disp('panels not ordered correctly');
    end

end

cx = sum(Cx);
cy = sum(Cy);

err_cx = abs(sum(fx_as(:))-cx);
err_cy = abs(sum(fy_as(:))-cy);

% check if total f_surface is equal to total Cp
if (err_cx >1e-8 || err_cy>1e-8 )
    disp(['integrated force is not conserved: ' num2str(err_cx) ' ' num2str(err_cy)]);
    
end

cl = -cx*sind(aoa) + cy*cosd(aoa)
cd = cx*cosd(aoa) + cy*sind(aoa)