function [Fx,Fy,dFx,dFy] = actuator_Maarten_force(V,t,options,getJacobian)
% force expression for 2D uniformly loaded actuator disk
% dFx, dFy are the Jacobians dFx/dV and dFy/dV

xin = options.grid.xin;
Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;

% Jacobian:
Nu  = options.grid.Nu;
Nv  = options.grid.Nv;
dFx = spalloc(Nu,Nu+Nv,0);
dFy = spalloc(Nv,Nu+Nv,0);    

indu = options.grid.indu;
indv = options.grid.indv;
uh   = V(indu);
vh   = V(indv);

order4 = options.discretization.order4;

% actuator disks with only forcing in x-direction, so acting on the
% (xin,yp) mesh
Fx = zeros(options.grid.Nu,1);
% vertical force is zero    
Fy = zeros(options.grid.Nv,1);

for k=1:length(options.force.x_c)
    
    Ct  = options.force.Ct(k);
    D   = options.force.D(k);
    x_c = options.force.x_c(k);
    y_c = options.force.y_c(k);
    
    % x-index of actuator disk
    x_ind = find(xin>=x_c,1);
    % y-index of actuator disk
%     y_ind = find(yp>=y_c,1);
    
    % disp(['position: ' num2str(xin(x_ind))])
    
    % get horizontal force
    % integrate in vertical direction
    y1Dy       = zeros(Nuy_in,1);
    for m=1:Nuy_in-1
        y1Dy(m) = integral(@(y) loading(y,y_c,D),options.grid.y(m),options.grid.y(m+1));
    end
    y1Dx        = zeros(Nux_in,1);
    y1Dx(x_ind) = 1;        
    
    if (order4==0)
        % based on U_inf
        % Fx = 0.5*(U_ref^2)*A*Ct
        % Ct = 4*a*(1-a)
%         Fx          = Fx - 0.5*(options.fluid.U1)^2*Ct*kron(y1Dy,y1Dx);
        % Jacobian will be zero in this case
        
        % based on local U
        % Fx = 0.5*(U_local^2)*A*Ct'
        % Ct' = Ct/(1-a)^2 = 4*a/(1-a)
        Fx          = Fx - 0.5*Ct*kron(y1Dy,y1Dx).*(uh.^2);
        if (getJacobian == 1)
           dFx      = dFx - [spdiags(0.5*Ct*kron(y1Dy,y1Dx).*(2.*uh),0,Nu,Nu) spalloc(Nu,Nv,0)]; 
        end

    elseif (order4==1)
        % fourth order:
        error('fourth order actuator disk not implemented');
%         
%         % small volumes
%         y1Dy1       = zeros(Nuy_in,1);
%         y1Dy1(yrange) = hy(yrange);
%         y1Dx1       = zeros(Nux_in,1);
%         y1Dx1(x_ind) = 1;
%         Fx1         = kron(y1Dy1,y1Dx1);
%         
%         % coarse volumes
%         y1Dy3         = zeros(Nuy_in,1);
%         y1Dy3(yrange(1)+1:yrange(end)-1) = 3;
%         y1Dy3(yrange(1)) = 2;
%         y1Dy3(yrange(end)) = 2;
%         y1Dy3(yrange(1)-1) = 1;
%         y1Dy3(yrange(end)+1) = 1;
%         
%         y1Dy3         = y1Dy3.*hy;
%         
%         y1Dx3_1       = zeros(Nux_in,1);
%         y1Dx3_1(x_ind-1:x_ind+1) = 1;
%         
%         Fx3           = kron(y1Dy3,y1Dx3_1);
%         
%         Fx            = -0.5*Ct*(alfa*Fx1 - Fx3);
    end
   
    
end

end

function f = loading(y,y_mid,D)
% get loading as a function of spatial coordinate
% used to integrate the loading over a cell
% y: vector
% y_mid, D: scalars

f = 1*(abs(y-y_mid)<=D/2);


end