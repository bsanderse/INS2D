% test integration of (1/r) over a wedge
% note that the clipping algorithm only works for CONVEX polygons (all interior
% angles<180 degrees), so general wedges are NOT supported (only if r1=0
% we have a convex polygon)

clear all
close all
% clc

tic

%% wedge
theta1 = 0;
theta2 = pi/2;

r1 = 0.;
r2 = 0.5;
x0 = 0.;
y0 = 0.;

% exact solutions ONLY hold for x0=y0=0 ?

% exact integral int int dr dtheta
A_ex = (r2-r1)*(theta2-theta1);
% exact integral int int -sin(theta) dr dtheta
B_ex = (r2-r1)*(cos(theta2)-cos(theta1));
% exact integral int int cos(theta) dr dtheta
C_ex = (r2-r1)*(sin(theta2)-sin(theta1));
% exact area int int r dr dtheta
S_ex = 0.5*(r2^2-r1^2)*(theta2-theta1);

%% construct a polygon as approximation to wedge
theta1 = linspace(theta1,theta2,10)';
theta2 = 0;%theta1;
x1 = x0+r2*cos(theta1);
x2 = x0+r1*cos(theta2);
x  = [x1;flipud(x2);x1(1)];
y1 = y0+r2*sin(theta1);
y2 = y0+r1*sin(theta2);
y  = [y1;flipud(y2);y1(1)];

subject = [x,y];


% plot(x,y,'rx-')
% axis([-1 1 -1 1]);
% grid
% axis equal

%% compute integral with polygon approximation
dy = y(2:end)-y(1:end-1);
dx = x(2:end)-x(1:end-1);
ds = sqrt(dx.^2+dy.^2);
nx = dy./ds;
ny = -dx./ds;
% int = sum(line_integralx(x(2:end),a,b)-line_integralx(x(1:end-1),a,b));
int = abs(sum(line_integral(x(2:end),y(2:end),nx,ny,ds,1)-line_integral(x(1:end-1),y(1:end-1),nx,ny,ds,1)) );


%% clip polygon onto mesh and integrate the different functions

% mesh
Nx = 10; Ny = 20;
xmesh = linspace(-2,2,Nx);
ymesh = linspace(-2,2,Ny);
% xmid  = 0.5*(xmesh(1:end-1)+xmesh(2:end));
% ymid  = 0.5*(ymesh(1:end-1)+ymesh(2:end));

Area  = zeros(Nx-1,Ny-1);
Area2 = zeros(Nx-1,Ny-1);
Area3 = zeros(Nx-1,Ny-1);
Area4 = zeros(Nx-1,Ny-1);
fun = @(x,y) 1;
fun2 = @(x,y) 1/sqrt(x^2+y^2);
fun3 = @(x,y) -y/(x^2+y^2);
for i=1:Nx-1
   
    for j=1:Ny-1

        % using the geom2D library
        clippedPolygon = clipPolygon(subject,[xmesh(i) xmesh(i+1) ymesh(j) ymesh(j+1)]);

        if (~isempty(clippedPolygon))

            % prepare for use of doubleintegral method (from file exchange)
            domain = struct('type','polygon','x',clippedPolygon(:,1)','y',clippedPolygon(:,2));
            param = struct('method','gauss','points',4);
            
            
            % area
%             Area(i,j) = abs(polygonArea(clippedPolygon));
%             Area(i,j) = polygonweightedArea(clippedPolygon,0);

            % weighted area: int int 1/r dx dy
            Area2(i,j) = polygonweightedArea(clippedPolygon,1);
%             Area2(i,j) = doubleintegral(fun2,domain,param);
        
            % weighted area: int int -y/r^2 dr dtheta
%             Area3(i,j) = polygonweightedArea(clippedPolygon,2);
%             Area3(i,j) = doubleintegral(fun3,domain,param,0,1);
%             pause

            
            % weighted area: int int x/r^2 dr dtheta
%             Area4(i,j) = polygonweightedArea(clippedPolygon,3);
        end
        
        % only for plotting:
%         clippingPolygon = [[xmesh(i);xmesh(i+1);xmesh(i+1);xmesh(i)],...
%                        [ymesh(j);ymesh(j);ymesh(j+1);ymesh(j+1)]];
%                    
        % plot
%         plot([subject(:,1);subject(1,1)],[subject(:,2);subject(1,2)],'blue')
%         hold on
%         patch(clippedPolygon(:,1),clippedPolygon(:,2),'rx-')
%         plot(clippingPolygon(:,1),clippingPolygon(:,2),'k')
%         pause;
        
    end
    
end

S_ex
A_ex
B_ex
C_ex

sum(sum(Area))
sum(sum(Area2))
sum(sum(Area3))
sum(sum(Area4))

toc

%% compute integral by integrating polygons on mesh

% Area = comp_Area(x',y',xmesh,ymesh,1);
% plot_staggered(xmesh,ymesh);
% hold on
% contour(xmid,ymid,Area);
% % hold on
% plot(x,y,'kx-')
% axis([0 1 0 1]);
% % grid
% axis equal
% 
% sum(sum(Area))