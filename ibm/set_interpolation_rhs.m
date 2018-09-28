function [rhs] = set_interpolation_rhs(xbi,ybi,xip,yip,xmesh,ymesh,interface_q,func,funcB)

% bi: boundary intercept
% ip: image point (in fluid)
% bp: bilinear points
% inside_vector(nq): vector with nq interface points, having coordinates
% index_vector(nq,1:2)

% func: function to be interpolated
% funcB: boundary value

% Nx   = length(xmesh);
% Ny   = length(ymesh);

nq   = length(xip);

rhs  = zeros(nq,4);

for q=1:nq

   
    indx1 = find(xmesh<xip(q),1,'last');
    indx2 = find(xmesh>=xip(q),1,'first');
    
    indy1 = find(ymesh<yip(q),1,'last');    
    indy2 = find(ymesh>=yip(q),1,'first');

    point = 1;
    
    % check for image points coinciding with mesh line
    
    % check for image points coinciding with geometry
    
    for i=indx1:indx2
        for j=indy1:indy2

            % find interface index q corresponding to point i,j
            q_indx = find(sum( abs(interface_q(:,1)-i) + abs(interface_q(:,2)-j),2) ==0);
%             keyboard
            if (q_indx>=1)
                % interface point, use value on body
                rhs(q,point) = funcB;                 
                
            else
                % fluid point
                rhs(q,point) = func(i,j);

            end
            
            point = point + 1;
        end
    end
    
%     if (point>5)
%         error('too many interpolation points used!');
%     end
 
%     keyboard        
%     
end