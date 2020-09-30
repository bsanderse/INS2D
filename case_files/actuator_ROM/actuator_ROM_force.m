function [Fx,Fy,dFx,dFy] = actuator_ROM_force(~,t,options,~)
% force expression for 2D uniformly loaded actuator disk

xin = options.grid.xin;
hy  = options.grid.hy;
Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;

% Jacobian:
Nu = options.grid.Nu;
Nv = options.grid.Nv;
dFx = spalloc(Nu,Nu+Nv,0);
dFy = spalloc(Nv,Nu+Nv,0);    


order4 = options.discretization.order4;

% actuator disk
Ct = options.force.Ct;
D  = options.force.D;
x_c = options.force.x_c;

% x-index of actuator disk
xmid = find(xin>=x_c,1);
 
% disp(['position: ' num2str(xin(xmid))])

Fy = zeros(Nvx_in,Nvy_in);
Fy = Fy(:);

if (rem(Nuy_in,2)==0)
    % Nuy_in even
    ymid   = Nuy_in/2;
    % number of volumes in disk (assuming locally uniform in y-dir)
    nd     = round(D/hy(ymid));
    yrange = ymid-floor(nd/2)+1:ymid+floor(nd/2);
else
    % Nuy_in odd
    ymid   = floor(Nuy_in/2)+1;
    % number of volumes in disk (assuming locally uniform in y-dir)
    nd     = D/hy(ymid);
    yrange = ymid-floor(nd/2):ymid+floor(nd/2);
end
    
    % construct a lift distribution resembling a real turbine blade
%     steepness   = 10;
%     if (rem(Nuy_in,2)==0)
%         coord       = yp(ymid+1:ymid+floor(nd/2));
%         force_distr = coord-D/2*(exp(steepness*coord)-1)/(exp(steepness*D/2) -1);
%         force_distr = [flipud(force_distr); force_distr];       
%     else
%         coord       = yp(ymid:ymid+floor(nd/2));
%         force_distr = coord-D/2*(exp(steepness*coord)-1)/(exp(steepness*D/2) -1);
%         force_distr = [flipud(force_distr(2:end)); force_distr];    
%     end
%     % normalize
%     const       = length(yrange)/sum(force_distr);

    
% ypos       = find(yp>=1,1);
% yrange     =  ypos - floor(nd/2)+1:ypos + floor(nd/2);

% full disc:
% yrange     = 1:Npy;

                    % uniform loading
% y1Dy(yrange)= const*force_distr.*hy(yrange);      % non-uniform loading

if (order4==0)
    y1Dy       = zeros(Nuy_in,1);
    y1Dy(yrange) = hy(yrange);          
    y1Dx       = zeros(Nux_in,1);
    y1Dx(xmid)  = 1;
    Fx          = -0.5*Ct*kron(y1Dy,y1Dx);
elseif (order4==1)
    % fourth order:

    % small volumes
    y1Dy1       = zeros(Nuy_in,1);
    y1Dy1(yrange) = hy(yrange);          
    y1Dx1       = zeros(Nux_in,1);
    y1Dx1(xmid) = 1;    
    Fx1         = kron(y1Dy1,y1Dx1);
    
    % coarse volumes
    y1Dy3         = zeros(Nuy_in,1);
    y1Dy3(yrange(1)+1:yrange(end)-1) = 3;
    y1Dy3(yrange(1)) = 2;
    y1Dy3(yrange(end)) = 2;
    y1Dy3(yrange(1)-1) = 1;
    y1Dy3(yrange(end)+1) = 1;
    
    y1Dy3         = y1Dy3.*hy;
    
    y1Dx3_1       = zeros(Nux_in,1);
    y1Dx3_1(xmid-1:xmid+1) = 1;    

    Fx3           = kron(y1Dy3,y1Dx3_1);
              
    Fx            = -0.5*Ct*(alfa*Fx1 - Fx3);
end

if (abs(sum(hy(yrange))-D)>1e-10)
    disp('total forcing not consistent, maybe the grid spacing cannot match the diameter');
    keyboard
end

Fx = Fx*(1 + sin(pi*t));

% Ct = 0.5;
% % second actuator
% y1Dx       = zeros(Nux_in,1);
% y1Dx(2*floor(Nx/3)) = 1;
% y1Dy       = zeros(Nuy_in,1);
% y1Dy(yrange+floor(nd/2))= hy(yrange);
% 
% Fx         = Fx -0.5*Ct*kron(y1Dy,y1Dx);

% 
% % third actuator
% y1Dx       = zeros(Nux_in,1);
% y1Dx(2*floor(Nx/5)) = 1;
% y1Dy       = zeros(Nuy_in,1);
% y1Dy(yrange-floor(Ny/4))= hy(yrange);
% 
% Fx         = Fx -0.5*Ct*kron(y1Dy,y1Dx);
% 
% % fourth actuator
% y1Dx       = zeros(Nux_in,1);
% y1Dx(3*floor(Nx/5)) = 1;
% y1Dy       = zeros(Nuy_in,1);
% y1Dy(yrange)= hy(yrange);
% 
% Fx         = Fx -0.5*Ct*kron(y1Dy,y1Dx);
% 
% % fifth actuator
% y1Dx       = zeros(Nux_in,1);
% y1Dx(4*floor(Nx/5)) = 1;
% y1Dy       = zeros(Nuy_in,1);
% y1Dy(yrange+floor(Ny/4))= hy(yrange);
% 
% Fx         = Fx -0.5*Ct*kron(y1Dy,y1Dx);
% 
% % sixth actuator
% y1Dx       = zeros(Nux_in,1);
% y1Dx(4*floor(Nx/5)) = 1;
% y1Dy       = zeros(Nuy_in,1);
% y1Dy(yrange-floor(Ny/4))= hy(yrange);
% 
% Fx         = Fx -0.5*Ct*kron(y1Dy,y1Dx);
    
end