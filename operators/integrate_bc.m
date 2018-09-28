% integrate BC
% used for RK methods

dudtLo_RK(:,i_RK) = dudtBC(x,y(1),t,Re); 
dudtUp_RK(:,i_RK) = dudtBC(x,y(end),t,Re);
dudtLe_RK(:,i_RK) = dudtBC(x(1),y,t,Re);
dudtRi_RK(:,i_RK) = dudtBC(x(end),y,t,Re);

uLo = uLo_n + dt*dudtLo_RK*A_RK(i_RK,:)';
uUp = uUp_n + dt*dudtUp_RK*A_RK(i_RK,:)';
uLe = uLe_n + dt*dudtLe_RK*A_RK(i_RK,:)';
uRi = uRi_n + dt*dudtRi_RK*A_RK(i_RK,:)';

dudtLo_RK_i(:,i_RK) = dudtBC(xin,y(1),t,Re); 
dudtUp_RK_i(:,i_RK) = dudtBC(xin,y(end),t,Re);
dudtLe_RK_i(:,i_RK) = dudtBC(x(1),yp,t,Re);
dudtRi_RK_i(:,i_RK) = dudtBC(x(end),yp,t,Re);

uLo_i = uLo_i_n + dt*dudtLo_RK_i*A_RK(i_RK,:)';
uUp_i = uUp_i_n + dt*dudtUp_RK_i*A_RK(i_RK,:)';
uLe_i = uLe_i_n + dt*dudtLe_RK_i*A_RK(i_RK,:)';
uRi_i = uRi_i_n + dt*dudtRi_RK_i*A_RK(i_RK,:)';

dvdtLo_RK(:,i_RK) = dvdtBC(x,y(1),t,Re); 
dvdtUp_RK(:,i_RK) = dvdtBC(x,y(end),t,Re);
dvdtLe_RK(:,i_RK) = dvdtBC(x(1),y,t,Re);
dvdtRi_RK(:,i_RK) = dvdtBC(x(end),y,t,Re);

vLo = vLo_n + dt*dvdtLo_RK*A_RK(i_RK,:)';
vUp = vUp_n + dt*dvdtUp_RK*A_RK(i_RK,:)';
vLe = vLe_n + dt*dvdtLe_RK*A_RK(i_RK,:)';
vRi = vRi_n + dt*dvdtRi_RK*A_RK(i_RK,:)';

dvdtLo_RK_i(:,i_RK) = dvdtBC(xp,y(1),t,Re); 
dvdtUp_RK_i(:,i_RK) = dvdtBC(xp,y(end),t,Re);
dvdtLe_RK_i(:,i_RK) = dvdtBC(x(1),yin,t,Re);
dvdtRi_RK_i(:,i_RK) = dvdtBC(x(end),yin,t,Re);

vLo_i = vLo_i_n + dt*dvdtLo_RK_i*A_RK(i_RK,:)';
vUp_i = vUp_i_n + dt*dvdtUp_RK_i*A_RK(i_RK,:)';
vLe_i = vLe_i_n + dt*dvdtLe_RK_i*A_RK(i_RK,:)';
vRi_i = vRi_i_n + dt*dvdtRi_RK_i*A_RK(i_RK,:)';