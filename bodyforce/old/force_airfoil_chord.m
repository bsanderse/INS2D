% prescribe Cl and Cd along airfoil chord

% staggered grid, so u- and v-volumes differ

% currently naca only
airfoil_type = '2412';
aoa          = 0;   % in angle of attack
Re_c         = 1e6; % Reynolds number or 0 for inviscid
xfoil        = 1;

% load airfoil properties, Cp and Cf from xfoil or own panel method
% nk panel endpoints, nk-1 panels
load_Cp_Cf;

% set coordinates of leading edge
% x_c       = L_x/3;
% y_c       = L_y/2;

x_k       = x_k + x_c;
y_k       = y_k + y_c;
x_cl      = x_cl + x_c;
y_cl      = y_cl + y_c;

% find leading edge and trailing edge, which determine chord
x_le      = x_k(ile);
y_le      = y_k(ile);
x_te      = x_k(1);
y_te      = y_k(1);
if (x_k(ile)==x_k(ile+1))
    y_le  = (y_k(ile)+y_k(ile+1))/2;
end
x_chord   = [x_le; x_te];
y_chord   = [y_le; y_te];

% find tangential and normal vectors of camberline
Sk        = sqrt((x_chord(2)-x_chord(1)).^2 + (y_chord(2)-y_chord(1)).^2);
vect(:,1) = (x_chord(2)-x_chord(1))/Sk;
vect(:,2) = (y_chord(2)-y_chord(1))/Sk;
vecn(:,1) = -vect(:,2);
vecn(:,2) = vect(:,1);


% find intersections of airfoil geometry with mesh.

%% u-volumes

[xi_u, yi_u, panel, param] = geometry_intersection(xp,y,x_chord,y_chord,Sk,0);

plot_staggered(xp,y)
hold on
plot(x_chord,y_chord,'rx-')
plot(xi_u,yi_u,'bo-')
plot(x_k,y_k,'gx-');

% for each finite volume intersected by a line segment find the integrated force
ni_u  = length(xi_u);
fx_as = zeros(Nux_in,Nuy_in);

% CpSx  = cx*vecn(:,1);
% CfSx  = cy*vect(:,1);

% a negative Cp is in the direction of the unit vector
% uniform distribution over entire chord:
Cx    = cx;

for i=1:ni_u-1
    k1 = panel(i);
    k2 = panel(i+1);
    
    % find indices of finite volume
    % points inside the FV
    
        xtemp = 0.5*(xi_u(i)+xi_u(i+1));
        ytemp = 0.5*(yi_u(i)+yi_u(i+1)); 

        % closest u-point=index FV
        [valx indx] = min(abs(xin-xtemp));
        [valy indy] = min(abs(yp-ytemp));
        
        fact  = 1;
        % check if minimum is unique
        % note: min() returns first value of minimum
        if (abs(xin(indx)-xtemp)==abs(xin(indx+1)-xtemp))
            indx = indx:indx+1;
            fact = fact+1;
        end
        if (abs(yp(indy)-ytemp)==abs(yp(indy+1)-ytemp))
            indy = indy:indy+1;  
            fact = fact+1;            
        end

    % find integrated force
    if (k2-k1==0)
      fx_as(indx,indy) = fx_as(indx,indy) + (1/fact)* ( Cx(k1)*(param(i+1)-param(i)) );
    elseif (k2-k1==1)
      fx_as(indx,indy) = fx_as(indx,indy) + (1/fact)* ( Cx(k1)*(1-param(i)) + Cx(k2)*param(i+1) );
    elseif (k2-k1>1)
      fx_as(indx,indy) = fx_as(indx,indy) + (1/fact)* ( Cx(k1)*(1-param(i)) + Cx(k2)*param(i+1) + ...
                         sum(Cx(k1+1:k2-1)) );
    else
      disp('panels not ordered correctly');    
    end
end

%% v-volumes
[xi_v, yi_v, panel, param] = geometry_intersection(x,yp,x_chord,y_chord,Sk,0);
% plot_staggered(x,yp)
% 
% hold on
% plot(x_cl,y_cl,'rx-')
% plot(xi_v,yi_v,'go-')

% for each finite volume intersected by a line segment find the integrated force
ni_v  = length(xi_v);
fy_as = zeros(Nvx_in,Nvy_in);

% CpSy  = C.*Sk.*vecn(:,2);
% CfSy  = Cf_cl.*Sk.*vect(:,2);

% a negative Cp is in the direction of the unit vector
Cy    = cy;

for i=1:ni_v-1
    k1 = panel(i);
    k2 = panel(i+1);
    
    % find indices of finite volume
    % points inside the FV
    
        xtemp = 0.5*(xi_v(i)+xi_v(i+1));
        ytemp = 0.5*(yi_v(i)+yi_v(i+1)); 

        % closest v-point=index FV
        [valx indx] = min(abs(xp-xtemp));
        [valy indy] = min(abs(yin-ytemp));

        fact  = 1;
        % check if minimum is unique; if not redistribute over neighbouring
        % volumes
        % note: min() returns first value of minimum
        if (abs(xp(indx)-xtemp)==abs(xp(indx+1)-xtemp))
            indx = indx:indx+1;
            fact = fact+1;
        end
        if (abs(yin(indy)-ytemp)==abs(yin(indy+1)-ytemp))
            indy = indy:indy+1;  
            fact = fact+1;            
        end
        
    % find integrated force
    if (k2-k1==0)
      fy_as(indx,indy) = fy_as(indx,indy) + (1/fact)* ( Cy(k1)*(param(i+1)-param(i)) );
    elseif (k2-k1==1)
      fy_as(indx,indy) = fy_as(indx,indy) + (1/fact)* ( Cy(k1)*(1-param(i)) + Cy(k2)*param(i+1) );
    elseif (k2-k1>1)
      fy_as(indx,indy) = fy_as(indx,indy) + (1/fact)* ( Cy(k1)*(1-param(i)) + Cy(k2)*param(i+1) + ...
                         sum(Cy(k1+1:k2-1)) );            
    else
      disp('panels not ordered correctly');
    end

end


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

keyboard;