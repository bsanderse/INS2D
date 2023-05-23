function PTybc = ybc_entries(ps,coords,t,options)

xs = coords(:,1);
ys = coords(:,2);

%%

xin = options.grid.xin;
% yin = options.grid.yin;
x   = options.grid.x;
% y   = options.grid.y;
% xp  = options.grid.xp;
yp  = options.grid.yp;

% Nxp = numel(xp);
Nx = numel(x);
Nxin = numel(xin);

Nyp = numel(yp);
% Ny = numel(y);
% Nyin = numel(yin);

% xs_u = xs(1:2*(Nx+Nxin+Nyp));
% xs_v = xs(1+2*(Nx+Nxin+Nyp):end);
% 
% ys_u = ys(1:2*(Nx+Nxin+Nyp));
% ys_v = ys(1+2*(Nx+Nxin+Nyp):end);
% 
% %%
% 
% us = uBC(xs_u,ys_u,t,options);
% vs = vBC(xs_v,ys_v,t,options);

global uBC vBC

% us = uBC(xs,ys,t,options); % uBC only accepts scalar input!
% vs = vBC(xs,ys,t,options); % vBC only accepts scalar input!

m = numel(ps);

us = zeros(m,1);
vs = zeros(m,1);

for l=1:m  % might be more expensive than get_bc_vector_ybc
    us(l) = uBC(xs(l),ys(l),t,options);
    vs(l) = vBC(xs(l),ys(l),t,options);
end


PTybc = zeros(m,1);
PTybc(ps<=2*(Nx+Nxin+Nyp)) = us(ps<=2*(Nx+Nxin+Nyp));
PTybc(ps >2*(Nx+Nxin+Nyp)) = vs(ps >2*(Nx+Nxin+Nyp));
