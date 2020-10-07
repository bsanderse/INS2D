%% THIS FILE IS NOW OBSOLETE AND REPLACED BY SET_BC_VECTORS

% make boundary values vectors if scalars
uLo = uLo.*ones(Nx+1,1);
uUp = uUp.*ones(Nx+1,1);
uLe = uLe.*ones(Ny+1,1);
uRi = uRi.*ones(Ny+1,1);

vLo = vLo.*ones(Nx+1,1);
vUp = vUp.*ones(Nx+1,1);
vLe = vLe.*ones(Ny+1,1);
vRi = vRi.*ones(Ny+1,1);

pLo = pLo.*ones(Nx+1,1);
pUp = pUp.*ones(Nx+1,1);
pLe = pLe.*ones(Ny+1,1);
pRi = pRi.*ones(Ny+1,1);

kLo = kLo.*ones(Nx+1,1);
kUp = kUp.*ones(Nx+1,1);
kLe = kLe.*ones(Ny+1,1);
kRi = kRi.*ones(Ny+1,1);

eLo = eLo.*ones(Nx+1,1);
eUp = eUp.*ones(Nx+1,1);
eLe = eLe.*ones(Ny+1,1);
eRi = eRi.*ones(Ny+1,1);


uLo_i      = uBC(options.grid.xin,y(1),t,options); 
uUp_i      = uBC(options.grid.xin,y(end),t,options); 
uLe_i      = uBC(x(1),options.grid.yp,t,options); 
uRi_i      = uBC(x(end),options.grid.yp,t,options); 

vLo_i      = vBC(options.grid.xp,y(1),t,options); 
vUp_i      = vBC(options.grid.xp,y(end),t,options);
vLe_i      = vBC(x(1),options.grid.yin,t,options);
vRi_i      = vBC(x(end),options.grid.yin,t,options);

% uLe_i  = interp1q(y,uLe,yp);
% uRi_i  = interp1q(y,uRi,yp);
% uLo_i  = interp1q(x,uLo,xin);
% uUp_i  = interp1q(x,uUp,xin);
% vLe_i  = interp1q(y,vLe,yin);
% vRi_i  = interp1q(y,vRi,yin);
% vLo_i  = interp1q(x,vLo,xp);
% vUp_i  = interp1q(x,vUp,xp);


% BC used when additional pressure equation is solved


dudtLe_i   = dudtBC(x(1),options.grid.yp,t,Re); 
dudtRi_i   = dudtBC(x(end),options.grid.yp,t,Re);

dvdtLo_i   = dvdtBC(options.grid.xp,y(1),t,Re); 
dvdtUp_i   = dvdtBC(options.grid.xp,y(end),t,Re);



% p, k and e are always needed at the xp, yp positions; so we interpolate
% it to these points beforehand
pLo = interp1q(x,pLo,xp);
pUp = interp1q(x,pUp,xp);
pLe = interp1q(y,pLe,yp);
pRi = interp1q(y,pRi,yp);

kLo = interp1q(x,kLo,xp);
kUp = interp1q(x,kUp,xp);
kLe = interp1q(y,kLe,yp);
kRi = interp1q(y,kRi,yp);

eLo = interp1q(x,eLo,xp);
eUp = interp1q(x,eUp,xp);
eLe = interp1q(y,eLe,yp);
eRi = interp1q(y,eRi,yp);