function [xi,yi,panel,param] = geometry_intersection(xmesh,ymesh,xk,yk,Sk,closed)
% computes intersection of a 2D Cartesian mesh with a 1D closed contour
% input: 
%  xmesh, ymesh are vectors containing coordinates of mesh lines
%  xk, yk are vectors containing coordinates of the contour
%  closed indicates if contour should be closed
% output:
%  xi, yi are vectors containing all intersection points, order by panel
%  number
%  panel is a vector containing the panel associated with an intersection
%  point
%  param is a vector containing the value of t for each intersection

% 
% imin = 1;
% imax = length(xmesh)-1;
% 
% jmin = 1;
% jmax = length(ymesh)-1;

% limit 'search' domain:
imin = max(find(xmesh>min(xk),1)-2,1);
imax = find(xmesh>max(xk),1)+1;
if (isempty(imax) || imax>length(xmesh)-1)
    imax = length(xmesh)-1;
end
jmin = max(find(ymesh>min(yk),1)-2,1);
jmax = find(ymesh>max(yk),1)+1;
if (isempty(jmax) || jmax>length(ymesh)-1)
    jmax = length(ymesh)-1;
end

% remove points outside domain?


nk   = length(xk)-1;

% loop over all sides of finite volumes


% initialize with number of intersections equal to number of panels
q     = 1;              % intersection counter
xi    = zeros(nk,1);    % intersection points
yi    = zeros(nk,1);
panel = zeros(nk,1);    % intersection panel
param = zeros(nk,1);    % parameter along panel
globalparam = zeros(nk,1);    % global parameter along contour


sumSk = cumsum(Sk);
sumSk = [0;sumSk];

% vertical faces
for i=imin:imax+1
    for j=jmin:jmax
        
            % s is the parameter ranging from 0 to 1 from a to b
            xa = xmesh(i);
            xb = xmesh(i);
            
            ya = ymesh(j);
            yb = ymesh(j+1);
            
            
            for k=1:nk

                % t is the parameter ranging from 0 to 1 from c to d
                xc = xk(k);
                yc = yk(k);
                
                xd = xk(k+1);
                yd = yk(k+1);

                D  = xa*(yd-yc) + xb*(yc-yd) + xd*(yb-ya) + xc*(ya-yb);
                if (D~=0)
                   s = (xa*(yd-yc) + xc*(ya-yd) + xd*(yc-ya))/D;
                   t = -(xa*(yc-yb) + xb*(ya-yc) + xc*(yb-ya))/D;
                   if (s>=0 && s<=1 && t>=0 && t<=1)

                       xi(q,1)    = xa + s*(xb-xa);
                       yi(q,1)    = ya + s*(yb-ya);
                       panel(q,1) = k;
                       param(q,1) = t;
                       globalparam(q,1) = sumSk(k) + t*Sk(k);
                       q     = q+1;
                   end

                end
                

            end

       
    end
end


% horizontal faces
for i=imin:imax
    for j=jmin:jmax+1             

            
            xa = xmesh(i);
            xb = xmesh(i+1);
            
            ya = ymesh(j);
            yb = ymesh(j);         
            
            for k=1:nk

                xc = xk(k);
                yc = yk(k);
                
                xd = xk(k+1);
                yd = yk(k+1);

                D  = xa*(yd-yc) + xb*(yc-yd) + xd*(yb-ya) + xc*(ya-yb);
                if (D~=0)
                   s = (xa*(yd-yc) + xc*(ya-yd) + xd*(yc-ya))/D;
                   t = -(xa*(yc-yb) + xb*(ya-yc) + xc*(yb-ya))/D;
                   if (s>=0 && s<=1 && t>=0 && t<=1)
                       
                       xi(q,1)    = xa + s*(xb-xa);
                       yi(q,1)    = ya + s*(yb-ya);
                       panel(q,1) = k;
                       param(q,1) = t; 
                       globalparam(q,1) = sumSk(k) + t*Sk(k);
                       q     = q+1;
                   end


                end
                

            end

       
    end
end
xi = xi(1:q-1,1);
yi = yi(1:q-1,1);

% find unique values
% [veci,indxi] = unique([xi yi],'rows');
% trick to use unique with round-off errors:
[veci,indxi] = unique(round([xi yi]*1e10),'rows');
veci  = 1e-10*veci;

xi    = veci(:,1);
yi    = veci(:,2);
panel = panel(indxi);
param = param(indxi);
globalparam = globalparam(indxi);

% sort according to globalparam 
[globalparam,indx] = sort(globalparam);
xi    = xi(indx);
yi    = yi(indx);
panel = panel(indx);
param = param(indx);

if (closed==1)
    % add last point to begin to get closed curve
    xi    = [xi(end); xi];
    yi    = [yi(end); yi];
    panel = [panel(end); panel];
    param = [param(end); param];
    eps   = 1e-12;
    if (abs(sum(diff(param)))>eps)
        disp(['parameter sum not equal to 0: ' num2str(sum(diff(param)))]);
    end
else
    eps   = 1e-12;    
    % add first and last point
    if (abs(xk(1)-xi(1))>eps || abs(yk(1)-yi(1))>eps )
        xi    = [xk(1); xi];
        yi    = [yk(1); yi];
        panel = [1; panel];
        param = [0; param];
    end
    if (abs(xk(end)-xi(end))>eps || abs(yk(end)-yi(end))>eps)
        xi    = [xi; xk(end)];
        yi    = [yi; yk(end)];        
        panel = [panel; nk];
        param = [param; 1];
    end
%     panel = [1; panel; nk];
%     param = [0; param; 1];
    if (abs(sum(diff(param))-1)>eps)
        disp(['parameter sum not equal to 1: ' num2str(sum(diff(param)))]);
    end  
end

