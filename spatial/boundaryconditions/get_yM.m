function yM = get_yM(t,options,yBC)
% construct boundary conditions


%% get settings from options structure

% steady
steady = options.case.steady;
% 4th order
order4 = options.discretization.order4;
% Reynolds number
Re = options.fluid.Re;
% boundary conditions
BC = options.BC;

global uBC vBC dudtBC dvdtBC;

% type of stress tensor
visc = options.case.visc;


if (order4 == 1)
    alfa   = options.discretization.alfa;
end

%% grid settings

% number of interior points and boundary points
Nux_in = options.grid.Nux_in;
Nvy_in = options.grid.Nvy_in;
Np     = options.grid.Np;
Npx    = options.grid.Npx;
Npy    = options.grid.Npy;


xin = options.grid.xin;
yin = options.grid.yin;
x   = options.grid.x;
y   = options.grid.y;
hx  = options.grid.hx;
hy  = options.grid.hy;
xp  = options.grid.xp;
yp  = options.grid.yp;


%% get BC values
xl = length(x);
yl = length(y);
xinl = length(xin);
yinl = length(yin);
xpl = length(xp);
ypl = length(yp);

uLo = yBC(1:xl); yBC = yBC(xl+1:end);
uUp = yBC(1:xl); yBC = yBC(xl+1:end);

uLo_i = yBC(1:xinl); yBC = yBC(xinl+1:end);
uUp_i = yBC(1:xinl); yBC = yBC(xinl+1:end);
uLe_i = yBC(1:ypl); yBC = yBC(ypl+1:end);
uRi_i = yBC(1:ypl); yBC = yBC(ypl+1:end);

vLe = yBC(1:yl); yBC = yBC(yl+1:end);
vRi = yBC(1:yl); yBC = yBC(y1+1:end);

vLo_i = yBC(1:xpl); yBC = yBC(xpl+1:end);
vUp_i = yBC(1:xpl); yBC = yBC(xpl+1:end);
vLe_i = yBC(1:yinl); yBC = yBC(yinl+1:end);
vRi_i = yBC(1:yinl); %yBC = yBC(yinl+1:end);

% uLo      = uBC(x,y(1),t,options);
% uUp      = uBC(x,y(end),t,options);
% % uLe      = uBC(x(1),y,t,options);
% % uRi      = uBC(x(end),y,t,options);
% 
% uLo_i    = uBC(xin,y(1),t,options);
% uUp_i    = uBC(xin,y(end),t,options);
% uLe_i    = uBC(x(1),yp,t,options);
% uRi_i    = uBC(x(end),yp,t,options);
% 
% % vLo      = vBC(x,y(1),t,options);
% % vUp      = vBC(x,y(end),t,options);
% vLe      = vBC(x(1),y,t,options);
% vRi      = vBC(x(end),y,t,options);
% 
% vLo_i    = vBC(xp,y(1),t,options);
% vUp_i    = vBC(xp,y(end),t,options);
% vLe_i    = vBC(x(1),yin,t,options);
% vRi_i    = vBC(x(end),yin,t,options);


pLe = options.BC.pLe;
pRi = options.BC.pRi;
pLo = options.BC.pLo;
pUp = options.BC.pUp;


%% boundary conditions for divergence

Mx_BC = options.discretization.Mx_BC;
My_BC = options.discretization.My_BC;


if (order4==1)
    Mx_BC3 = options.discretization.Mx_BC3;
    My_BC3 = options.discretization.My_BC3;
end


% Mx
ybc    = kron(uLe_i,Mx_BC.ybc1) + kron(uRi_i,Mx_BC.ybc2);
yMx    = Mx_BC.Bbc*ybc;
if (order4==1)
    ybc3 = kron(uLe_i,Mx_BC3.ybc1) + kron(uRi_i,Mx_BC3.ybc2);
    yMx3 = Mx_BC3.Bbc*ybc3;
    yMx  = alfa*yMx - yMx3;
end

% My
ybc    = kron(My_BC.ybc1,vLo_i) + kron(My_BC.ybc2,vUp_i);
yMy    = My_BC.Bbc*ybc;
if (order4==1)
    ybc3 = kron(My_BC3.ybc1,vLo_i) + kron(My_BC3.ybc2,vUp_i);
    yMy3 = My_BC3.Bbc*ybc3;
    yMy  = alfa*yMy - yMy3;
end

yM     = yMx + yMy;
% options.discretization.yM  = yM;

