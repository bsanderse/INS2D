function [Area] = comp_Area(xx,yy,x,y,weight)
%--------------------------------------------------------------------------
%  Computes Area's of intersection of contour (xx,yy) with mesh (x,y)
%--------------------------------------------------------------------------
% if weight==1 then we compute int 1/r dx dy instead of int dx dy
if (nargin==4)
    weight=0;
end   

Nx = length(x);
Ny = length(y);
Nelem = length(xx)-1;

% Close contour
xx = [xx(end), xx, xx(1)];
yy = [yy(end), yy, yy(1)];

% Confine contour
yL = y(1:end-1)';
yR = y(2:end)';

x2 = zeros((Nx-1)*(Ny-1), 2*(2*Nelem+2)+2);
y2 = zeros((Nx-1)*(Ny-1), 2*(2*Nelem+2)+2);
for i=1:Nx-1
%   x-dir
%   fprintf('i=%i\n', i)
  xL = x(i);
  xR = x(i+1);

  [x1,y1] = confine_x(xx,yy,xL,xR);
%   figure(1)
%   plot(xx,yy, x1,y1)
%   fprintf('Press a key...\n')
%   keyboard
  
% Confine contour y-dir
  [ydum,xdum] = confine_x(y1,x1,yL,yR);

  idx = i:Nx-1:i+(Ny-2)*(Nx-1);
  x2(idx,:) = xdum;
  y2(idx,:) = ydum;
end
Area = reshape(comp_area_single(x2,y2,weight), Nx-1, Ny-1)';
% keyboard

function [xdum, ydum] = confine_x(xx,yy, xL, xR)
%--------------------------------------------------------------------------
%  
%--------------------------------------------------------------------------
Nelem = length(xx(:))-1;
Nxy   = length(xL(:));

xL = repmat(reshape(xL, Nxy, 1), 1, Nelem);
xR = repmat(reshape(xR, Nxy, 1), 1, Nelem);

x1 = repmat(reshape(xx(1:end-1), 1, Nelem), Nxy, 1);
x2 = repmat(reshape(xx(2:end), 1, Nelem), Nxy, 1);

y1 = repmat(reshape(yy(1:end-1), 1, Nelem), Nxy, 1);
y2 = repmat(reshape(yy(2:end), 1, Nelem), Nxy, 1);

fdum = (y2-y1).*(x2-x1)./((x2-x1).^2+1E-16);

xi1 = (min(max(x1,xL),xR)-x1);
xi2 = (min(max(x2,xL),xR)-x1);

% reorder
%idx = find(abs(xi2)<abs(xi1));
%dum = xi1(idx);
%xi1(idx) = xi2(idx);
%xi2(idx) = dum;

%idx = find(xi1==xi2);
%xi1(idx) = NaN;
%xi2(idx) = NaN;

xp1=x1+xi1;
xp2=x1+xi2;

yp1=y1+xi1.*fdum;
yp2=y1+xi2.*fdum;

% xdum = [xp1(1:end); xp2(1:end)];
% ydum = [yp1(1:end); yp2(1:end)];

xdum = [reshape(xp1',1,Nxy*Nelem); reshape(xp2',1,Nxy*Nelem)];
ydum = [reshape(yp1',1,Nxy*Nelem); reshape(yp2',1,Nxy*Nelem)];

 
% fprintf('(%5.2f %5.2f) -> (%5.2f %5.2f)\n', [xp1; yp1; xp2; yp2])

xdum = reshape(xdum, 2*(Nelem), Nxy)';
ydum = reshape(ydum, 2*(Nelem), Nxy)';

function [area] = comp_area_single(xx,yy,weight)
%-----------------------------------------------------
%  computes area enclosed by (xx,yy)
%-----------------------------------------------------
% if weight = 0 we compute the area, if weight=1 we compute int 1/r dx dy

if (weight==0)
    xL = xx(:,1:end-1);
    yL = yy(:,1:end-1);

    xR = xx(:,2:end);
    yR = yy(:,2:end);

    tx = xR-xL;
    ty = yR-yL;

    area = abs(sum(0.5*(xL.*ty-yL.*tx),2));
elseif (weight==1)
  
%     [M,N]=size(xx);
%     area = zeros(M,1);
%     for i=1:M
%     dy = yy(i,2:end)-yy(i,1:end-1);
%     dx = xx(i,2:end)-xx(i,1:end-1);
%     a  = dy./dx;
%     b  = yy(i,1:end-1)-a.*xx(i,1:end-1);     
%     
%     area(i)=sum(line_integralx(xx(i,2:end),a,b)-line_integralx(xx(i,1:end-1),a,b));
%     
%     end
    
    % vector form:
%     dy = yy(:,2:end)-yy(:,1:end-1);
%     dx = xx(:,2:end)-xx(:,1:end-1);
%     a  = dy./dx;
%     b  = yy(:,1:end-1)-a.*xx(:,1:end-1);     
    dy = yy(:,2:end)-yy(:,1:end-1);
    dx = xx(:,2:end)-xx(:,1:end-1);
    ds = sqrt(dx.^2+dy.^2);
    nx = dy./ds;
    ny = -dx./ds;    
%     nx(abs(ds)<eps)=0; ny(abs(ds)<eps)=0;
%     keyboard;
    area = abs(sum(line_integralx(xx(:,2:end),yy(:,2:end),nx,ny,ds) ...
                  -line_integralx(xx(:,1:end-1),yy(:,1:end-1),nx,ny,ds),2));
    
end