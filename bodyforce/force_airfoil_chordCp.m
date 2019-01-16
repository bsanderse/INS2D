% prescribe Cp and Cf along airfoil camber

% staggered grid, so u- and v-volumes differ

% currently naca only
airfoil_type = '2412';
aoa          = 4;   % in angle of attack
Re_c         = 5e2; % Reynolds number or 0 for inviscid
xfoil        = 1;
rotated      = 1;

% load airfoil properties, Cp and Cf from xfoil or own panel method
% nk panel endpoints, nk-1 panels
% result: 
% xfoil: Cp_ep and Cf_ep; pressure and skin friction at panel end points (x_k,y_k)
%        interpolate to: Cp and Cf at 'collocation' points (xcol,ycol) (midpoints)
% own method: Cp and Cf at xcol,ycol

load_Cp_Cf;

% the camberline from Xfoil is not so smooth
% redefine the camberline as follows:
% -take same number of points in x-dir on upper and lower side
% -interpolate y
% -average upper and lower side

nk_cl     = 100; 
idx       = (1:nk_cl+1)';
x_ku      = x_k(1:ile);   % upper side, from TE to LE
x_kl      = x_k(ile:end); % lower side, from LE to TE
y_ku      = y_k(1:ile);
y_kl      = y_k(ile:end);

% resample new x-coordinates (panel end points),
% equal on both lower and upper side (from TE
% to LE). this also equals the camberline distribution
% here we assume that y and Cp are functions of x (i.e. each x has a unique
% y and Cp)
x_rs_l    = x_k(ile)+(x_k(1)-x_k(ile))*sin(0.5*pi*(idx-1)*(1/nk_cl)).^2;
x_rs_u    = flipud(x_rs_l); % from TE to LE
x_rs      = [x_rs_u; x_rs_l(2:end)];
% midpoints:
% x_rs_col  = 0.5*(x_rs_l(1:end-1)+x_rs_l(2:end));
% y-values
y_rs_u    = interp1(x_ku,y_ku,x_rs_u,'pchip','extrap');
y_rs_l    = interp1(x_kl,y_kl,x_rs_l,'pchip','extrap');
y_rs      = [y_rs_u; y_rs_l(2:end)];


% midpoints (alternatively one could directly interpolate to x_rs_col
% y_rs_col_u = 0.5*(y_rs_u(1:end-1)+y_rs_u(2:end));
% y_rs_col_l = 0.5*(y_rs_l(1:end-1)+y_rs_l(2:end));
% y_rs_col_u = interp1(x_ku,y_ku,x_rs_col,'pchip','extrap');
% y_rs_col_l = interp1(x_kl,y_kl,x_rs_col,'pchip','extrap');

% get normal and tangential vectors of resampled panels
% normal=outward normal; tangential vector points in direction of TE on
% both upper and lower side
[x_rs_col,y_rs_col,S_rs,vect_rs,vecn_rs] = panel_properties(x_rs,y_rs,1,1);
x_rs_col_u = x_rs_col(1:nk_cl);
x_rs_col_l = x_rs_col(nk_cl+1:end);
y_rs_col_u = y_rs_col(1:nk_cl);
y_rs_col_l = y_rs_col(nk_cl+1:end);
% visualize with:
% plotnormal(x_rs,y_rs,vect_rs,vecn_rs)

% Cp and Cf values at end points
% Cp_rs_u   = interp1(x_ku,Cp_ep(1:ile),x_rs,'pchip','extrap');
% Cp_rs_l   = interp1(x_kl,Cp_ep(ile:end),x_rs,'pchip','extrap');
% Cf_rs_u   = interp1(x_ku,Cf_ep(1:ile),x_rs,'pchip','extrap');
% Cf_rs_l   = interp1(x_kl,Cf_ep(ile:end),x_rs,'pchip','extrap');
% Cp and Cf values at midpoints of resampled panels
Cp_rs_u   = interp1(x_ku,Cp_ep(1:ile),x_rs_col_u,'pchip','extrap');
Cp_rs_l   = interp1(x_kl,Cp_ep(ile:end),x_rs_col_l,'pchip','extrap');
Cp_rs     = [Cp_rs_u; Cp_rs_l];
Cf_rs_u   = interp1(x_ku,Cf_ep(1:ile),x_rs_col_u,'pchip','extrap');
Cf_rs_l   = interp1(x_kl,Cf_ep(ile:end),x_rs_col_l,'pchip','extrap');
Cf_rs     = [Cf_rs_u; Cf_rs_l];

% compare resampled total force with resampled pressure distribution to
% original total force
% Cp_rs = Cp_rs*sum(Cp.*Sk)/sum(Cp_rs.*S_rs);
% Cf_rs = Cf_rs*sum(Cf.*Sk)/sum(Cf_rs.*S_rs);
disp('original total force')
sum(-Cp.*Sk.*vecn(:,1) + Cf.*Sk.*vect(:,1)) % x-dir
sum(-Cp.*Sk.*vecn(:,2) + Cf.*Sk.*vect(:,2)) % y-dir
disp('resampled total force')
sum(-Cp_rs.*S_rs.*vecn_rs(:,1) + Cf_rs.*S_rs.*vect_rs(:,1))
sum(-Cp_rs.*S_rs.*vecn_rs(:,2) + Cf_rs.*S_rs.*vect_rs(:,2))

% the coordinates of the chord line; from TE to LE
x_cl      = x_rs_u;
y_cl      = (y_k(1)-y_k(ile))/(x_k(1)-x_k(ile))*x_cl + y_k(ile);
% y_cl      = [0.5*(y_rs_u(1)+ y_rs_l(end)); 0.5*(y_rs_u(end)+ y_rs_l(1))];
% the point at the leading edge can be troublesome
% we obtain it here as an extrapolation from the interior points
% y_cl(end)   = interp1(x_cl(1:end-1),y_cl(1:end-1),x_cl(end),'pchip','extrap');
[x_cl_col,y_cl_col,S_cl,vect_cl,vecn_cl] = panel_properties(x_cl,y_cl,0,1);

% plotnormal(x_cl,y_cl,vect_cl,vecn_cl);

if (rotated==1)
    [x_cl,y_cl] = rotate_body(x_cl,y_cl,aoa);
end
        
x_cl = x_cl + x_c;
y_cl = y_cl + y_c;


% transfer forces to chordline
% x-component
Cx_cl     = -Cp_rs(1:nk_cl).*S_rs(1:nk_cl).*vecn_rs(1:nk_cl,1) + ...
            -Cp_rs(end:-1:nk_cl+1).*S_rs(end:-1:nk_cl+1).*vecn_rs(end:-1:nk_cl+1,1) + ...
            Cf_rs(1:nk_cl).*S_rs(1:nk_cl).*vect_rs(1:nk_cl,1) + ...
            Cf_rs(end:-1:nk_cl+1).*S_rs(end:-1:nk_cl+1).*vect_rs(end:-1:nk_cl+1,1);
% y-component
Cy_cl     = -Cp_rs(1:nk_cl).*S_rs(1:nk_cl).*vecn_rs(1:nk_cl,2) + ...
            -Cp_rs(end:-1:nk_cl+1).*S_rs(end:-1:nk_cl+1).*vecn_rs(end:-1:nk_cl+1,2) + ...
            Cf_rs(1:nk_cl).*S_rs(1:nk_cl).*vect_rs(1:nk_cl,2) + ...
            Cf_rs(end:-1:nk_cl+1).*S_rs(end:-1:nk_cl+1).*vect_rs(end:-1:nk_cl+1,2);

if (rotated==1) % lift and drag are aligned with x-axis
    Cx_cl_rot = Cx_cl*cosd(aoa) + Cy_cl*sind(aoa);
    Cy_cl_rot = -Cx_cl*sind(aoa) + Cy_cl*cosd(aoa);
    Cx_cl = Cx_cl_rot;
    Cy_cl = Cy_cl_rot;
end
    
% now find intersections of airfoil geometry with mesh.


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

