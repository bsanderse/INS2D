
function plot_staggered(x,y,z)
% plot velocity and pressure points for staggered grid
% argument z is optional
% x,y and z are position of velocity points ('grid lines'), 
% i.e. u-velocities lie on the x=const. lines, v-velocities lie on the
% y=const lines, w-velocities lie on the z=const lines

% note: Nx = no. of finite volumes in x-dir + 1
%       Ny = no. of finite volumes in y-dir + 1
%       Nz = no. of finite volumes in z-dir + 1

 if (nargin==2)
    
    Nx = length(x);
    Ny = length(y);

    figure
    hold on

    % grid lines
    for i=1:Nx
       plot([x(i) x(i)],[y(1) y(Ny)],'k-');     
    end

    for j=1:Ny
       plot([x(1) x(Nx)],[y(j) y(j)],'k-'); 
    end

    % u-velocity points
%     for i=1:Nx
%         for j=1:Ny-1
%             plot(x(i),(y(j+1)+y(j))/2,'bx');
%         end
%     end
% 
%     % v-velocity points
%     for j=1:Ny
%         for i=1:Nx-1
%             plot((x(i+1)+x(i))/2,y(j),'rx');
%         end
%     end
% 
%     % pressure points
%     for i=1:Nx-1
%        for j=1:Ny-1
%            plot((x(i+1)+x(i))/2,(y(j+1)+y(j))/2,'go');
%        end
%     end

% axis equal

 elseif (nargin==3)
   
    Nx = length(x);
    Ny = length(y);
    Nz = length(z);
    
    figure
    hold on

    for k=1:Nz
        for i=1:Nx
            plot3([x(i) x(i)],[y(1) y(Ny)],[z(k) z(k)],'k-');     
        end            
    end

    for k=1:Nz
        for j=1:Ny
            plot3([x(1) x(Nx)],[y(j) y(j)],[z(k) z(k)],'k-'); 
        end
    end   
    
    for j=1:Ny
        for i=1:Nx
            plot3([x(i) x(i)],[y(j) y(j)],[z(1) z(Nz)],'k-'); 
        end
    end
    
    view(3);
    axis([x(1) x(end) y(1) y(end) z(1) z(end)]);
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
 end
 
end
