imin=1;
imax=Nx;

jmin=1;
jmax=Ny;

% loop over all sides of finite volumes
% initialize

q = 1;                  % intersection counter
xi = zeros(nk,1);       % intersection points
yi = zeros(nk,1);
panel = zeros(nk,1);    % intersection panel
param = zeros(nk,1);    % parameter along panel

% vertical faces

for i=imin:imax+1
    for j=jmin:jmax
        
            % s is the parameter ranging from 0 to 1 from a to b
            xa = x(i);
            xb = x(i);
            
            ya = y(j);
            yb = y(j+1);
            
            
            for k=1:nk

                % t is the parameter ranging from 0 to 1 from c to d
                xc = x_k(k);
                yc = y_k(k);
                
                xd = x_k(k+1);
                yd = y_k(k+1);

                D  = xa*(yd-yc) + xb*(yc-yd) + xd*(yb-ya) + xc*(ya-yb);
                if (D~=0)
                   s = (xa*(yd-yc) + xc*(ya-yd) + xd*(yc-ya))/D;
                   t = -(xa*(yc-yb) + xb*(ya-yc) + xc*(yb-ya))/D;
                   if (s>=0 && s<=1 && t>=0 && t<=1)

                       xi(q) = xa + s*(xb-xa);
                       yi(q) = ya + s*(yb-ya);
                       panel(q) = k;
                       param(q) = t;
                     
                       q     = q+1;
                   end

                end
                

            end

       
    end
end

% horizontal faces
for i=imin:imax
    for j=jmin:jmax+1
        
       
            xa = x(i);
            xb = x(i+1);
            
            ya = y(j);
            yb = y(j);
            
            xc = x_k(1:nk-1);
            yc = y_k(1:nk-1);
            xd = x_k(2:nk);
            yd = y_k(2:nk);
            
            D  = xa*(yd-yc) + xb*(yc-yd) + xd*(yb-ya) + xc*(ya-yb);
            % note: D can be zero, leading to Inf values in s and t
                          
            s = (xa.*(yd-yc) + xc.*(ya-yd) + xd.*(yc-ya))./D;
            t = -(xa.*(yc-yb) + xb.*(ya-yc) + xc.*(yb-ya))./D;
            
            test  = s>=0 & s<=1 & t>=0 & t<=1;
            % get out zeros with unique
            xi    = unique(test.*(xa + s.*(xb-xa)));
            yi    = unique(text.*(ya + s.*(yb-ya)));
            
            xi    = xi(2:end);
            yi    = yi(2:end);
            
            panel(q) = k;
            param(q) = t; 

%             
%             if (s>=0 && s<=1 && t>=0 && t<=1)
%             
%                    xi(q) = xa + s*(xb-xa);
%                    yi(q) = ya + s*(yb-ya);
%                    panel(q) = k;
%                    param(q) = t; 
% 
%                    q     = q+1;
% 
%             end            
            
%             for k=1:nk
% 
%                 xc = x_k(k);
%                 yc = y_k(k);
%                 
%                 xd = x_k(k+1);
%                 yd = y_k(k+1);
% 
%                 D  = xa*(yd-yc) + xb*(yc-yd) + xd*(yb-ya) + xc*(ya-yb);
%                 if (D~=0)
%                    s = (xa*(yd-yc) + xc*(ya-yd) + xd*(yc-ya))/D;
%                    t = -(xa*(yc-yb) + xb*(ya-yc) + xc*(yb-ya))/D;
%                    if (s>=0 && s<=1 && t>=0 && t<=1)
%                        
%                        xi(q) = xa + s*(xb-xa);
%                        yi(q) = ya + s*(yb-ya);
%                        panel(q) = k;
%                        param(q) = t; 
% 
%                        q     = q+1;
%                    end
% 
% 
%                 end
%                 
% 
%             end

       
    end
end
xi = xi(1:q-1);
yi = yi(1:q-1);

% find unique values
[veci,indxi] = unique([xi yi],'rows');
xi   = veci(:,1);
yi   = veci(:,2);
panel = panel(indxi);
param = param(indxi);

% sort according to panel number
[panel,indx] = sort(panel);
xi   = xi(indx);
yi   = yi(indx);
param= param(indx);

% add first point to end to get closed curve
xi = [xi(end); xi];
yi = [yi(end); yi];
panel = [1; panel];
param = [0; param];