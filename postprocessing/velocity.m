% velocity values at pressure points

[up,vp,qp] = get_velocity(V,t,options);

% BFS:
% l = [0.05 0.1 0.15 0.2 0.4 0.6 0.8 1 1.2 1.4];
% list = 0.1:0.05:1.3;
list = 25;

figure(10)
% pcolor(xp,yp,qp')
contour(xp,yp,qp',list)
% surf(xp,yp,qp');
% contour(xp,yp,up',25);
% shading interp
axis([x1 x2 y1 y2]);
 axis equal
% axis tight
% axis square
colorbar

% print(gcf,'-depsc',['results/nonaligned_convection/surf_' num2str(Nx)]);

xlabel('$x$','Interpreter','Latex');
ylabel('$y$','Interpreter','Latex');


%% get wake profiles
% u = reshape(uh,Nux_in,Nuy_in);
% v = reshape(vh,Nvx_in,Nvy_in);
% uwake1 = interp2(xin',yp,u',x_c+0.5,yp);
% uwake2 = interp2(xin',yp,u',x_c+5,yp);
% vwake1 = interp2(xp',yin,v',x_c+0.5,yin);
% vwake2 = interp2(xp',yin,v',x_c+5,yin);
% 
% 
% figure(2)
% plot(uwake1,yp,'rx-')
% hold on
% plot(uwake2,yp,'ro-')
% figure(3)
% plot(vwake1,yin,'bx-')
% hold on
% plot(vwake2,yin,'bo-')


% figure(2)
% contour(xp,yp,qp',25)
% axis equal
% axis tight
% axis([x1 x2 y1 y2])
% colorbar
% print(gcf,'-depsc',['results/nonaligned_convection/contour_' num2str(Nx)]);


% actuator disk
% hold on
% h1=plot([x(xmid+1) x(xmid+1)],[yp(yrange(1)) yp(yrange(end))],'k-');
% set(h1,'LineWidth',1.5);
% h2=plot([x(2*xmid+1) x(2*xmid+1)],[yp(yrange(1)+floor(nd/2)) yp(yrange(end)+floor(nd/2))],'k-');
% set(h2,'LineWidth',1.5);

% airfoil
% hold on
% plot(x_k,y_k,'bx-');

% hold off;
