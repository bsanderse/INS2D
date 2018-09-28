% prescribe Cp and Cf along airfoil surface

% staggered grid, so u- and v-volumes differ

% currently naca only
airfoil_type = '2412';
aoa          = 4;   % in angle of attack
Re_c         = 5e2; % Reynolds number or 0 for inviscid
xfoil        = 1;
rotated      = 1;

% load airfoil properties, Cp and Cf from xfoil or own panel method
% load airfoil properties, Cp and Cf from xfoil or own panel method
% nk panel endpoints, nk-1 panels
% result: 
% xfoil: Cp_ep and Cf_ep; pressure and skin friction at panel end points (x_k,y_k)
%        interpolate to: Cp and Cf at 'collocation' points (xcol,ycol) (midpoints)
% own method: Cp and Cf at xcol,ycol
load_Cp_Cf;
    
if (rotated==1)
    [x_k,y_k] = rotate_body(x_k,y_k,aoa);
end
x_k       = x_k + x_c;
y_k       = y_k + y_c;

% plot tangential and normal vectors:
% figure;plotnormal(x_k,y_k,vect,vecn);

CpSx  = Cp.*Sk.*vecn(:,1);
CfSx  = Cf.*Sk.*vect(:,1);
CpSy  = Cp.*Sk.*vecn(:,2);
CfSy  = Cf.*Sk.*vect(:,2);

% a negative Cp is in the direction of the unit vector
Cx    = -CpSx+CfSx;
Cy    = -CpSy+CfSy;

if (rotated==1) % lift and drag are aligned with x-axis
    Cx_rot = Cx*cosd(aoa) + Cy*sind(aoa);
    Cy_rot = -Cx*sind(aoa) + Cy*cosd(aoa);
    Cx = Cx_rot;    
    Cy = Cy_rot;    
end
% find intersections of airfoil geometry with mesh.

%% u-volumes

[xi_u, yi_u, panel, param] = geometry_intersection(xp,y,x_k,y_k,Sk,0);
% plot_staggered(xp,y)
% hold on
% plot(x_k,y_k,'rx-')
% plot(xi_u,yi_u,'bo-')

% for each finite volume intersected by a line segment find the integrated force
fx_as = regularize_force(xin,yp,xi_u,yi_u,panel,param,Cx);


%% v-volumes
[xi_v, yi_v, panel, param] = geometry_intersection(x,yp,x_k,y_k,Sk,0);
% plot_staggered(x,yp)
% hold on
% plot(x_k,y_k,'rx-')
% plot(xi_v,yi_v,'go-')

% for each finite volume intersected by a line segment find the integrated force
fy_as = regularize_force(xp,yin,xi_v,yi_v,panel,param,Cy);


%%

Fx = -0.5*fx_as(:);
Fy = -0.5*fy_as(:);

%% checks
cx = sum(Cx);
cy = sum(Cy);

err_cx = abs(sum(fx_as(:))-cx);
err_cy = abs(sum(fy_as(:))-cy);

% check if total f_surface is equal to total Cp
if (err_cx >1e-8 || err_cy>1e-8 )
    disp(['integrated force is not conserved: ' num2str(err_cx) ' ' num2str(err_cy)]);
    
end
