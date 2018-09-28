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
xADu = [xmu1;xmu2]; yADu = [ymu1;ymu2];
xADv = [xmv1;xmv2]; yADv = [ymv1;ymv2];

% xAD = [xm1+hx(1)/2;xm2-hx(end)/2]; yAD = [ym1-hy(end)/2;ym2+hy(1)/2];
SADu = sqrt( diff(xADu)^2 + diff(yADu)^2);
SADv = sqrt( diff(xADv)^2 + diff(yADv)^2);

% force in x- and y-dir
Cx  = Ct*cosd(alfa)*SADu;
Cy  = Ct*sind(alfa)*SADv;

%% intersections with u-volumes
% boundaries of u-volumes given by xp,y
[xi_u,yi_u,panel,param] = geometry_intersection(xp,y,xADu,yADu,SADu,0);
% figure(1)
% plot_staggered(xp,y)
% hold on
% plot(xADu,yADu,'rx-')
% plot(xi_u,yi_u,'bo-')

fx_as = regularize_force(xin,yp,xi_u,yi_u,panel,param,Cx);
if (abs(sum(fx_as(:))-Cx)>1e-10)
    warning('force not conserved');
end

%% intersections with v-volumes
% boundaries of v-volumes given by x,yp
[xi_v,yi_v,panel,param] = geometry_intersection(x,yp,xADv,yADv,SADv,0);
% figure(2)
% plot_staggered(x,yp)
% hold on
% plot(xADv,yADv,'rx-')
% plot(xi_v,yi_v,'bo-')

fy_as = regularize_force(xp,yin,xi_v,yi_v,panel,param,Cy);
if (abs(sum(fy_as(:))-Cy)>1e-10)
    warning('force not conserved');
end

%% intersection with p-volumes for correction terms
% boundaries of u-volumes
% loop over all horizontal pressure faces:
AD = [xmv1 ymv1 xmu2 ymu2];
% figure
% plot_staggered(x,y)
% plot([AD(1) AD(3)],[AD(2) AD(4)])
% hold on

% fracu = zeros(Npx,Npy);
% fracv = zeros(Npx,Npy);
% 
% for ix=1:Npx
%     for jx=1:Npy
%    
%         faceu = [x(ix) yp(jx) x(ix+1) yp(jx)];
%         
%         pointu = intersectEdges(faceu,AD);
%         
%         if (~isnan(pointu))
% %             plot(pointu(1),pointu(2),'rx')
%             fracu(ix,jx) = sign(0.5*(x(ix)+x(ix+1))-pointu(1))*min(abs(pointu(1)-x(ix)),abs(pointu(1)-x(ix+1)))/deltax;
%         end
%         
%         facev = [xp(ix) y(jx) xp(ix) y(jx+1)];
%         
%         pointv = intersectEdges(facev,AD);
%         
%         if (~isnan(pointv))
% %             plot(pointv(1),pointv(2),'kx')
%             fracv(ix,jx) = sign(0.5*(y(jx)+y(jx+1))-pointv(2))*min(abs(pointv(2)-y(jx)),abs(pointv(2)-y(jx+1)))/deltay;            
%         end        
%         
%     end
% end
% Fx_corr = 0.5*Ct*Gx*fracu(:);
% Fy_corr = 0.5*Ct*Gy*fracv(:);

%%
Fx_corr = zeros(Nux_in,Nuy_in);
Fy_corr = zeros(Nvx_in,Nvy_in);

slope = (AD(4)-AD(2))/(AD(3)-AD(1));

plot_staggered(x,y)
plot([AD(1) AD(3)],[AD(2) AD(4)])
hold on
for ix=1:Nux_in-1
    for jx=1:Nuy_in
   
        face_mid = [xp(ix) yp(jx); xp(ix+1) yp(jx)];
        right_of_AD = +(face_mid(:,1)*slope+AD(4)-face_mid(:,2)<0); % logical->double

        if (length(unique(right_of_AD))>1)
        plot([face_mid(1,1) face_mid(2,1)],[face_mid(1,2) face_mid(2,2)],'rx-')
        Fx_corr(ix,jx) = hy(jx);
        
        end  
        
    end
end
for ix=1:Nux_in
    for jx=1:Nuy_in-1
   
        face_mid = [xp(ix) yp(jx); xp(ix) yp(jx+1)];
        right_of_AD = +(face_mid(:,1)*slope+AD(4)-face_mid(:,2)<0); % logical->double

        if (length(unique(right_of_AD))>1)
        plot([face_mid(1,1) face_mid(2,1)],[face_mid(1,2) face_mid(2,2)],'gx-')
        Fy_corr(ix,jx) = hx(ix);
        
        end  
        
    end
end
Fx = -0.5*Ct*Fx_corr(:);
Fy = -0.5*Ct*Fy_corr(:);

%%
% Fx = -0.5*fx_as(:) + Fx_corr(:);
% Fy = -0.5*fy_as(:) + Fy_corr(:);


%% finite difference type correction
AD = [xmv1 ymv1 xmu2 ymu2];
weight = zeros(Npx,Npy);
slope = (AD(4)-AD(2))/(AD(3)-AD(1));
% y= slope*x+AD(4)
% plot([AD(1) AD(3)],[AD(2) AD(4)])
% hold on
diagA = reshape(diag(M*spdiags(1./Om,0,Nu+Nv,Nu+Nv)*G),Npx,Npy);

for ix=1:Npx
    for jx=1:Npy
   
        stencil = [xp(ix) yp(jx)];
        if (ix>1)
            stencil = [stencil; xp(ix-1) yp(jx)];
        end
        if (ix<Npx)
            stencil = [stencil; xp(ix+1) yp(jx)];
        end
        if (jx>1)
            stencil = [stencil; xp(ix) yp(jx-1)];
        end
        if (jx<Npy)
            stencil = [stencil; xp(ix) yp(jx+1)];
        end

        right_of_AD = +(stencil(:,1)*slope+AD(4)-stencil(:,2)<0); % logical->double
        if (length(unique(right_of_AD))>1)
           right_of_AD(1) = right_of_AD(1)*diagA(ix,jx);
           weight(ix,jx) = sum(right_of_AD);
%            plot(xp(ix),yp(jx),'ro')
%            keyboard
        end
        
%         dist = distancePointEdge(stencil,AD)
        
    end
end
     
weight = -weight*0.5*Ct;


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



