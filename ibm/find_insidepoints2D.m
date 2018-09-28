function [interface_q, inside, inside_excl_interface] = find_insidepoints2D(x,y,xk,yk,cp)
% interface_q is a VECTOR with interface points
% inside is a MATRIX with the points in the body set to 1, 
% inside_excl_interface is the same as inside, but EXCLUDES the interface
% points

% x,y: coordinates of finite volume center
% xk,yk: body geometry
% cp: control point used for ray tracing

if (nargin<4)
    disp('at least 4 input arguments required');
end   
if (nargin<5)
    cp = [x(1) y(1)];
end


eps  = 1e-12;
Nx   = length(x);
Ny   = length(y);

% factor 2 is for safety
imin = find(x<min(xk),1,'last')-2;
imax = find(x>max(xk),1,'first')+2; 
jmin = find(y<min(yk),1,'last')-2; 
jmax = find(y>max(yk),1,'first')+2;

if (isempty(imin))
    imin=1;
end
if (isempty(imax))
    imax=Nx;
end
if (isempty(jmin))
    jmin=1;
end
if (isempty(jmax))
    jmax=Ny;
end

nk   = length(xk)-1;

%% find intersections


xa = cp(1);
ya = cp(2);
% 
inside = zeros(Nx,Ny);

for i=imin:imax
    for j=jmin:jmax
        
            % s is the parameter ranging from 0 to 1 from a to b
            xb = x(i);
            yb = y(j);
            
            % intersection counter            
            q     = 0;
            for k=1:nk

                % t is the parameter ranging from 0 to 1 from c to d
                xc = xk(k);
                yc = yk(k);
                
                xd = xk(k+1);
                yd = yk(k+1);

                D  = xa*(yd-yc) + xb*(yc-yd) + xd*(yb-ya) + xc*(ya-yb);
                % D=0 means parallel lines
                if (abs(D)>eps)
                   s = (xa*(yd-yc) + xc*(ya-yd) + xd*(yc-ya))/D;
                   t = -(xa*(yc-yb) + xb*(ya-yc) + xc*(yb-ya))/D;
                   
                   % 'normal' intersection
                   if (s>eps && s<1-eps && t>eps && t<1-eps)       
                       
                       q = q+1;

                   end
                   % vertex intersection                   
                   if ( (s>eps && s<1-eps) && abs(t)<=eps) 
                       % we can skip t=1, it will lead to double counting
                       % check if ray is not 'touching' body
%                        disp('vertex intersection')
                       % normal to the ray:
                       n_R   = [-(yb-ya);xb-xa];
                       if ( (k == nk-1 && t==1) || k==nk || ...
                            (k == 2 && t==0) || k==1 )
                           t_k1 = [xk(2)-xc;yk(2)-yc];                               
                           t_k2 = [xk(nk-1)-xc;yk(nk-1)-yc];
                       else
                           t_k1 = [xd-xc;yd-yc];                               
                           t_k2 = [xk(k-1)-xc;yk(k-1)-yc];                           
                       end
                               
                       sign1 = sign(n_R'*t_k1);
                       sign2 = sign(n_R'*t_k2);
                       
                       if (sign1~=sign2)
                           q = q+1;
                       end
                   end
                   if (abs(s-1)<=eps && t>eps && t<1-eps)
                       disp('mesh point lying exactly on geometry')

                   end
                   if (abs(s-1)<=eps && (abs(t)<=eps || abs(t-1)<=eps))
                       disp('mesh point lying on vertex')

                   end                   
                   if (abs(s)<=eps && (t>=-eps && t<=1+eps))
                       disp('choose control point further from geometry')
                   end
                else
                    % segment parallel to ray
                end
%                 keyboard


            end
            inside(i,j) = mod(q,2);

            

       
    end
end


% inside points are now 1, outside 0

% figure
% pcolor(x,y,inside');
% axis equal
% keyboard

%% find interface points

inside_excl_interface = inside;
% outerface = inside;

% counter for interface points, which are a dimension lower than the mesh
q = 1;
for i=imin:imax
    for j=jmin:jmax

        ind = [i-1 j; i+1,j; i,j-1; i,j+1];

        if (i==1)
            if (j==1)
                ind([1 3],:) = [];
            elseif (jmin==Ny)
                ind([1 4],:) = [];
            else
                ind(1,:) = [];
            end
        elseif (i==Nx)
            if (j==1)
                ind([2 3],:) = [];
            elseif (j==Ny)
                ind([2 4],:) = [];
            else
                ind(2,:) = [];
            end
        elseif (i>1 && i<Nx)
            if (j==1)
                ind(3,:) = [];
            elseif (j==Ny)
                ind(4,:) = [];
            end
        end
        ind = sub2ind([Nx Ny],ind(:,1),ind(:,2));
%         sum_surr           = sum(inside(i-1,j)+inside(i+1,j) +...
%                                  inside(i,j-1)+inside(i,j+1));
        sum_surr = sum(inside(ind));
        max_sum  = length(inside(ind));
        
        if (inside(i,j)==1 && sum_surr>0 && sum_surr<max_sum)
            % store i,j value 
            interface_q(q,1:2)  = [i,j];
            % set inside_excl_interface to 0 for this point
            inside_excl_interface(i,j) = 0;
            q = q + 1;
        end
%         if (sum_surr>0)
%             outerface(i,j) = 1;
%         end

    end
end