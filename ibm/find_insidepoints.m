function inside_matrix = find_insidepoints(x,y,xb,yb,D,z0)

% xmin= 0.2;
% xmax= 1.8;
% ymin= 0.1;
% ymax= 0.9;
% 
% imin= max(find(x>=xmin,1)-1,1);
% imax= min(find(x>=xmax,1),Nx);
% jmin= max(find(y>=ymin,1)-1,1);
% jmax= min(find(y>=ymax,1),Ny);

Nx   = length(x);
Ny   = length(y);

imin = 2;
imax = Nx-1;
jmin = 2;
jmax = Ny-1;

n    = length(xb);

%% midpoints and normals of body 'panels'
% xm  = (xb(1:end-1)+xb(2:end))/2;
% ym  = (yb(1:end-1)+yb(2:end))/2;


% determine surface normals
% equation for a triangle is: ax+by+cz=d, or in vector form:
% n*x = d, where x=(x,y,z) and n=(a,b,c)

% normal vector for a triangle:
A = [xb(1:end-1); yb(1:end-1); z0*ones(1,n-1)];
B = [xb(2:end); yb(2:end); z0*ones(1,n-1)];
C = B + [zeros(1,n-1);zeros(1,n-1); 1*ones(1,n-1)]; % a vector perpendicular to this line

normals = cross(B-A,C-A); % not normalized, so not a unit vector

d = dot(normals,A);       % which is the same as dot(normals,B)

%% find intersections

% control point
Dn = repmat(D,1,n-1);
% 
inside = zeros(Nx,Ny);

% for i=imin:imax
%     
%     for j=jmin:jmax
%         
% 
%         % query point
%         Q  = [x(i); y(j); z0];
%         Qn = repmat(Q,1,n-1);
% 
%         % ray
%         rn= Qn-Dn;
% 
%         % intersection of ray with a line segment?
%         % intersection point between ray and plane containing line segment;
%     
%         % parameter along ray:
%         % vectorized:
%         s  = (d - dot(Dn,normals)) ./ dot(rn,normals);
%         % coordinates of intersection point:
%         un = Dn + [rn(1,:).*s; rn(2,:).*s; rn(3,:).*s]; 
%             
%         Area  = cross(B-A,C-A);
%         Area  = sqrt(Area(1,:).^2 + Area(2,:).^2 + Area(3,:).^2);
%         Area2 = cross(C-B,un-B);
%         Area2 = sqrt(Area2(1,:).^2 + Area2(2,:).^2 + Area2(3,:).^2);
%         Area3 = cross(A-C,un-C);
%         Area3 = sqrt(Area3(1,:).^2 + Area3(2,:).^2 + Area3(3,:).^2);
% 
%         n_int = sum( (abs(s)<1).*((Area2+Area3-Area)<eps) );
% 
%         % even amount of intersections: outside; odd amount: inside
%         inside(i,j) = mod(n_int,2);
% 
%     end
%     
% end


for i=imin:imax
    
    for j=jmin:jmax
        
        
    % query point
    Q  = [x(i); y(j); z0];

    % ray
    r = Q-D;

    
    % intersection of ray with a line segment?
    % intersection point between ray and plane containing line segment
    n_int = 0;
    
    for t=1:n-1
        % parameter along ray:
        s = (d(t) - dot(D,normals(:,t))) / dot(r,normals(:,t));
        
        if (abs(s)<1)
            % coordinates of intersection point:
            u  = D + r*s; 
            
            % check if intersection point is 'real', i.e. it lies on the line segment
            % determine areas of ABU, BCU and CAU
            Area = norm(cross(B(:,t)-A(:,t),C(:,t)-A(:,t))); % length of segment
    %         Area1 = cross(B(:,t)-A(:,t),u-A(:,t)); % should be zero for 2D
            Area2 = norm(cross(C(:,t)-B(:,t),u-B(:,t)));
            Area3 = norm(cross(A(:,t)-C(:,t),u-C(:,t)));


            if (abs(Area2+Area3-Area)<1e-14)
                % check if the intersection point is not a point where the
                % ray is 'touching' the geometry
                u1   = D + r*(s+eps); 
                u2   = D + r*(s-eps);
                vol1 = det([A(:,t)' 1; B(:,t)' 1; C(:,t)' 1; u1' 1]);
                vol2 = det([A(:,t)' 1; B(:,t)' 1; C(:,t)' 1; u2' 1]);
                if (sign(vol1)~=sign(vol2))
                    n_int = n_int+1;
%                     keyboard
                end
            end
        end
        
    end
    inside(i,j) = mod(n_int,2);

   end
    
end

% inside points are now 1, outside 0

figure
pcolor(x,y,inside');
axis equal

%% find interface points
sum_surr      = zeros(Nx,Ny);
inside_matrix = zeros(Nx,Ny);

% counter for interface points, which are a dimension lower than the mesh
q = 1;
for i=imin:imax
    for j=jmin:jmax
%         sum_surr(i,j)      = (inside(i,j)==0)*sum(sum(inside(i-1:i+1,j-1:j+1)));
%         inside_matrix(i,j) = (sum_surr(i,j)>0 && inside(i,j)==0);
        qi(q,1) = i;
        qj(q,2) = j;
        sum_surr(i,j)      = (inside(i,j)==1)*sum(inside(i-1,j)+inside(i+1,j) +...
                                                  inside(i,j-1)+inside(i,j+1));
        inside_matrix(i,j) = sum_surr(i,j)>0 && sum_surr(i,j)<4;
        q = q + 1;
    end
end

% inside_matrix = inside_matrix + inside;