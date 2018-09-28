% prescribe Cl and Cd along airfoil chord

% staggered grid, so u- and v-volumes differ

% currently naca only
airfoil_type = '2412';
aoa          = 4;   % in angle of attack
Re_c         = 5e2; % Reynolds number or 0 for inviscid
xfoil        = 1;
rotated      = 1;

Cl = 0.2246;
Cd = 0.17267;

% load airfoil properties, Cp and Cf from xfoil or own panel method
% nk panel endpoints, nk-1 panels
% result: 
% xfoil: Cp_ep and Cf_ep; pressure and skin friction at panel end points (x_k,y_k)
%        interpolate to: Cp and Cf at 'collocation' points (xcol,ycol) (midpoints)
% own method: Cp and Cf at xcol,ycol

% load_Cp_Cf;

% the coordinates of the chord line; from TE to LE
% obtain from airfoil coordinates
% x_cl      = [x_k(1); x_k(ile)];
% y_cl      = (y_k(1)-y_k(ile))/(x_k(1)-x_k(ile))*x_cl + y_k(ile);

x_cl = [1;0];
y_cl = [0;0];

% only one panel, so this is really simple:
[x_cl_col,y_cl_col,S_cl,vect_cl,vecn_cl] = panel_properties(x_cl,y_cl,0,1);

% plotnormal(x_cl,y_cl,vect_cl,vecn_cl);

% plot(x_k,y_k,'rx-')
% hold on
% plot(x_cl,y_cl,'bo-')
% keyboard;

% find intersections of airfoil geometry with mesh.
if (rotated==1)
    [x_cl,y_cl] = rotate_body(x_cl,y_cl,aoa);
end

x_cl = x_cl + x_c;
y_cl = y_cl + y_c;

if (rotated==1) % lift and drag are aligned with x-axis and freestream
    Cx_cl = Cd;
    Cy_cl = Cl;
else
    CX    = Cd*cosd(aoa) - Cl*sind(aoa);
    CY    = Cd*sind(aoa) + Cl*cosd(aoa);
    Cx_cl = CX;
    Cy_cl = CY;    
end
    
%% u-volumes
% boundaries of u-volumes are given by (xp, y)
% get the sequence of points that make up the intersection of the contour
% with mesh, and for each point the corresponding panel number of the
% contour, and the parameter value along it at which intersection occurs
[xi_u, yi_u, panel, param] = geometry_intersection(xp,y,x_cl,y_cl,S_cl,0);
% plot_staggered(xp,y)
% hold on
% plot(x_cl,y_cl,'rx-')
% plot(xi_u,yi_u,'bo-')

fx_as = regularize_force(xin,yp,xi_u,yi_u,panel,param,Cx_cl);

%% v-volumes
% boundaries of u-volumes are given by (x, yp)
[xi_v, yi_v, panel, param] = geometry_intersection(x,yp,x_cl,y_cl,S_cl,0);
% plot_staggered(x,yp)
% 
% hold on
% plot(x_cl,y_cl,'rx-')
% plot(xi_v,yi_v,'go-')

fy_as = regularize_force(xp,yin,xi_v,yi_v,panel,param,Cy_cl);


%% 
% force on flow. factor (1/2) takes into account non-dimensionalization
Fx = -0.5*fx_as(:);
Fy = -0.5*fy_as(:);

%% checks

err_cx = abs(sum(fx_as(:))-sum(Cx_cl));
err_cy = abs(sum(fy_as(:))-sum(Cy_cl));

% check if total f_surface is equal to total Cp
if (err_cx >1e-8 || err_cy>1e-8 )
    disp(['integrated force is not conserved: ' num2str(err_cx) ' ' num2str(err_cy)]);
    
end

% as a check with Xfoil:
% note that errors have been introduced in the interpolation from endpoints
% to midpoints
% cl = -cx*sind(aoa) + cy*cosd(aoa);
% cd = cx*cosd(aoa) + cy*sind(aoa);

