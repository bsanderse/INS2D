% test clipping

% subject = [[-2; 0.5; 0.5; -2; -2],[-2; -2; 2; 2; -2]];
subject = [...
   0.841470984807897  -0.500000000000000;...
   0.997494986604055  -0.500000000000000;...
   0.997494986604055   0.500000000000000;...
   0.841470984807897   0.500000000000000;...
   0.841470984807897  -0.500000000000000];

% mesh = [[0; 1; 1; 0; 0],[0; 0; 1; 1; 0];

clippedPolygon = clipPolygon(subject,[0.98 1.02 -0.48 -0.44])


% clippingPolygon = [[0; 1; 1; 0; 0],[-1; -1; 0; 0; -1]];

figure
plot(subject(:,1),subject(:,2),'bx-')
hold on
patch(clippedPolygon(:,1),clippedPolygon(:,2),'rx-')
% plot(clippingPolygon(:,1),clippingPolygon(:,2),'k')
