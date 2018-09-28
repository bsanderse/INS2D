function [bilin_mat_inv, bilin_points] = set_interpolation_matrix(xbi,ybi,xip,yip,xmesh,ymesh,interface_q,neumann,xk,yk,panel_bi)

% bi: boundary intercept
% ip: image point (in fluid)
% bp: bilinear points
% inside_vector(nq): vector with nq interface points, having coordinates
% index_vector(nq,1:2)

% neumann: if 1 then use neumann boundary condition. default 0.

if (nargin<8)
    neumann = 0;
end

if (neumann == 1)
    [xcol,ycol,Sk,vect,vecn] = panel_properties(xk,yk,1,0);
end

Nx = length(xmesh);
Ny = length(ymesh);
nq = length(xip);

bilin_mat_inv = zeros(nq,4,4); % inverse of interpolation matrix for each q
bilin_points  = zeros(nq,4); % index values of interpolation points for each q

for q=1:nq

    x_bp  = zeros(4,1);
    y_bp  = zeros(4,1);

    bilin_mat = zeros(4,4);
    
    indx1 = find(xmesh<xip(q),1,'last');
    indx2 = find(xmesh>=xip(q),1,'first');
    
    indy1 = find(ymesh<yip(q),1,'last');    
    indy2 = find(ymesh>=yip(q),1,'first');

    point = 1;
    
    % check for image points coinciding with mesh line
    
    % check for image points coinciding with geometry
    
    for i=indx1:indx2
        for j=indy1:indy2

            % find if there is an interface index q corresponding to point
            % i,j
            q_indx = find(sum( abs(interface_q(:,1)-i) + abs(interface_q(:,2)-j),2) ==0);
            if (q_indx>=1)

                % interface point, use point on body
                x_bp(point) = xbi(q_indx);     
                y_bp(point) = ybi(q_indx);                
                               
                if (neumann ==0 )
                    bilin_mat(point,1:4) = ...
                        [1 xbi(q_indx) ybi(q_indx) xbi(q_indx)*ybi(q_indx)];
                else
                    vecnx = vecn(panel_bi(q_indx),1);
                    vecny = vecn(panel_bi(q_indx),2);
                    bilin_mat(point,1:4) = ...
                        [0 vecnx vecny xbi(q_indx)*vecny+ybi(q_indx)*vecnx];
                end
                    

            else
                % fluid point
                x_bp(point) = xmesh(i);
                y_bp(point) = ymesh(j);               
    
                bilin_points(q,point) = sub2ind([Nx,Ny],i,j);
                
                bilin_mat(point,1:4) = ...
                    [1 xmesh(i) ymesh(j) xmesh(i)*ymesh(j)];
            end
            point = point + 1;
        end
    end
    
%     if (point>5)
%         error('too many interpolation points used!');
%     end
    
    % inverse of 4x4 matrix
    bilin_mat_inv(q,1:4,1:4) = inv(bilin_mat);
    

%     plot(xip(q),yip(q),'rx');
%     hold on
%     plot(x_bp,y_bp,'mo');
%     keyboard;
%     
            
    
end