% actuator disk
Ct = 0.5;   % thrust coefficient
% D  = 2/cosd(30);     % diameter
%  

alfa = 30; % angle wrt vertical in degrees

% actuator coordinates
% xADu = [x_c - (D/2)*sind(alfa); x_c + (D/2)*sind(alfa)];
% yADu = [y_c + (D/2)*cosd(alfa); y_c - (D/2)*cosd(alfa)];
% xADv = xADu; yADv = yADu;
% SADu = sqrt( diff(xADu)^2 + diff(yADu)^2);
% SADv = SADu;

% coordinates 1D behind turbine
xm = -0.5; %D*cosd(alfa);
ym = -1+0.5*sqrt(3); %D*sind(alfa); 
% xm   = x1+cosd(alfa)^2;
% ym   = y1+cosd(alfa)*sind(alfa);
xmu1 = max(xp(1),xm-(y2-ym)*tand(alfa)); xmu2 = min(xp(end),xm+(ym-y1)*tand(alfa));
ymu1 = min(y2,ym-(xp(1)-xm)/tand(alfa)); ymu2 = max(y1,ym-(xp(end)-xm)/tand(alfa));
xmv1 = max(x1,xm-(yp(end)-ym)*tand(alfa)); xmv2 = min(x2,xm+(ym-yp(1))*tand(alfa));
ymv1 = min(yp(end),ym-(x1-xm)/tand(alfa)); ymv2 = max(yp(1),ym-(x2-xm)/tand(alfa));
% 
% 
Nk = Nx;
xADu = linspace(xmu1,xmu2,Nk); yADu = linspace(ymu1,ymu2,Nk);
xADu_mid = (xADu(1:end-1)+xADu(2:end))/2;
yADu_mid = (yADu(1:end-1)+yADu(2:end))/2;

xADv = linspace(xmv1,xmv2,Nk); yADv = linspace(ymv1,ymv2,Nk);
xADv_mid = (xADv(1:end-1)+xADv(2:end))/2;
yADv_mid = (yADv(1:end-1)+yADv(2:end))/2;

SADu = sqrt( diff(xADu).^2 + diff(yADu).^2);
SADv = sqrt( diff(xADv).^2 + diff(yADv).^2);

% force in x- and y-dir
fkx  = Ct*cosd(alfa)*SADu;
fky  = Ct*sind(alfa)*SADv;


%% 
% Dirac regularization
% naive implementation (slow)
% assuming uniform grid.
e   = 3*hx(1); 
Fx  = zeros(Nux_in,Nuy_in);
for i=1:Nux_in
    for j=1:Nuy_in
        
        for k=1:Nk-1
            
            rx = xin(i)-xADu_mid(k);
            ry = yp(j)-yADu_mid(k);
            dx = (hx(1)/(2*e))*(cos(pi*(rx/e))+1)*(abs(rx/e)<=1);
            dy = (hy(1)/(2*e))*(cos(pi*(ry/e))+1)*(abs(ry/e)<=1);
%             dx=(hx(1)/e)*(1/sqrt(pi))*exp(-(rx/e)^2);
%             dy=(hy(1)/e)*(1/sqrt(pi))*exp(-(ry/e)^2);

            Fx(i,j) = Fx(i,j)+dx*dy*fkx(k);
        end
        
    end
end

%%
e   = 3*hy(1); 
Fy  = zeros(Nvx_in,Nvy_in);
for i=1:Nvx_in
    for j=1:Nvy_in
        
        for k=1:Nk-1
            
            rx = xp(i)-xADv_mid(k);
            ry = yin(j)-yADv_mid(k);
            dx = (hx(1)/(2*e))*(cos(pi*(rx/e))+1)*(abs(rx/e)<=1);
            dy = (hy(1)/(2*e))*(cos(pi*(ry/e))+1)*(abs(ry/e)<=1);
%             dx=(hx(1)/e)*(1/sqrt(pi))*exp(-(rx/e)^2);
%             dy=(hy(1)/e)*(1/sqrt(pi))*exp(-(ry/e)^2);

            Fy(i,j) = Fy(i,j)+dx*dy*fky(k);
        end
        
    end
end

%%
Fx = -0.5*Fx(:); 
Fy = -0.5*Fy(:);


%% plotting stuff
% Fx2 = zeros(Nux_in+2,Nuy_in+2); % include boundary cells
% Fx2(2:Nux_in+1,2:Nuy_in+1) = reshape(Fx,Nux_in,Nuy_in);
% % Fx2 = Fx;
% Fy2 = Fy;
% % Fy2(abs(Fy2)<eps)=NaN;
% 
% figure
% pcolor([x(1);xin;xin(end)+hx(end)]-0.5*deltax,[y(1);yp;yp(end)+hy(end)]-0.5*deltay,Fx2')
% colormap('gray')
% hold on
% % plot(xAD,yAD,'rx-','Linewidth',2)
% plot(xi_u,yi_u,'rx-','Linewidth',2)
% axis([x1 x2 y1 y2])
% box on
% 
% figure
% pcolor(xp-0.5*deltax,yin-0.5*deltay,reshape(Fy2,Nvx_in,Nvy_in)')
% colormap('gray')
% hold on
% % plot(xAD,yAD,'rx-','Linewidth',2)
% plot(xi_v,yi_v,'rx-','Linewidth',2)



