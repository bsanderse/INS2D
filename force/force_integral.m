function Area = force_integral(subject,xmesh,ymesh,bnd,type)
% subject: polygon to be clipped
% xmesh, ymesh: coordinates of mesh lines
% bnd: boundaries to limit computation
% type: 
% 0=area, 
% 1=weighted with 1/r
% 2=weighted with -y/r^2
% 3=weighted with x/r^2

Nx = length(xmesh);
Ny = length(ymesh);

Area = zeros(Nx-1,Ny-1);
imin = bnd(1);
imax = bnd(2);
jmin = bnd(3);
jmax = bnd(4);

for i=imin:imax
   
    for j=jmin:jmax

        % using the geom2D library
        % find intersectional polygon of finite volume and geometry
        clippedPolygon = clipPolygon(subject,[xmesh(i) xmesh(i+1) ymesh(j) ymesh(j+1)]);

        
        % weighted area: int int 1/r dr dtheta
        if (~isempty(clippedPolygon))

%         clippingPolygon = [[xmesh(i);xmesh(i+1);xmesh(i+1);xmesh(i)],...
%                            [ymesh(j);ymesh(j);ymesh(j+1);ymesh(j+1)]];
%                    
%         plot([subject(:,1);subject(1,1)],[subject(:,2);subject(1,2)],'blue')
%         hold on
%         patch(clippedPolygon(:,1),clippedPolygon(:,2),'rx-')
%         plot(clippingPolygon(:,1),clippingPolygon(:,2),'k')

          if (type==0)
            Area(i,j) = polygonArea(clippedPolygon);
          else
            Area(i,j) = polygonweightedArea(clippedPolygon,type);
          end
       
        end
        
    end
    
end

%         % area
%         if (~isempty(clippedPolygon))
%         Area(i,j) = abs(polygonArea(clippedPolygon));
%         end
%         
%         % weighted area: int int 1/r dr dtheta
%         if (~isempty(clippedPolygon))
%         Area2(i,j) = polygonweightedArea(clippedPolygon,1);
%         end
%         
%         % weighted area: int int -y/r^2 dr dtheta
%         if (~isempty(clippedPolygon))
%         Area3(i,j) = polygonweightedArea(clippedPolygon,2);
%         end
%         
%         % weighted area: int int x/r^2 dr dtheta
%         if (~isempty(clippedPolygon))
%         Area4(i,j) = polygonweightedArea(clippedPolygon,3);
%         end
%         
