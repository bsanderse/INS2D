function coords = ybc_coords(ps,t,options)
% [xs_u,ys_u,xs_v,ys_v]

% m = numel(ps);

% xs_u = zeros(m,1);
% ys_u = zeros(m,1);
% xs_v = zeros(m,1);
% ys_v = zeros(m,1);

xin = options.grid.xin;
yin = options.grid.yin;
x   = options.grid.x;
y   = options.grid.y;
xp  = options.grid.xp;
yp  = options.grid.yp;

Nxp = numel(xp);
Nx = numel(x);
Nxin = numel(xin);

Nyp = numel(yp);
Ny = numel(y);
Nyin = numel(yin);

Xs_u = [x; x; xin; xin; x(1)*ones(Nyp,1); x(end)*ones(Nyp,1)];
Ys_u = [y(1)*ones(Nx,1); y(end)*ones(Nx,1); ...
        y(1)*ones(Nxin,1); y(end)*ones(Nxin,1); yp; yp];
Xs_v = [x(1)*ones(Ny,1); x(end)*ones(Ny,1); xp; xp; ...
        x(1)*ones(Nyin,1); x(end)*ones(Nyin,1)];
Ys_v = [y; y; y(1)*ones(Nxp,1); y(end)*ones(Nxp,1); yin; yin];

Xs = [Xs_u; Xs_v];
Ys = [Ys_u; Ys_v];

% xs = Xs(ps);
% ys = Ys(ps);
% 
% xs_u = xs(1:2*(Nx+Nxin+Nyp));
% xs_v = xs(1+2*(Nx+Nxin+Nyp):end);
% 
% ys_u = ys(1:2*(Nx+Nxin+Nyp));
% ys_v = ys(1+2*(Nx+Nxin+Nyp):end);

coords = [Xs(ps) Ys(ps)];