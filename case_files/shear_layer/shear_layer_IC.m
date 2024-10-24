function [u,v,p,options] = shear_layer_IC(t,options)
% initial velocity field

yu = options.grid.yu;
xv = options.grid.xv;
Npx = options.grid.Npx;
Npy = options.grid.Npy;

d   = pi/15;
e   = 0.05;
% note: we add 1 to u in order to make global momentum conservation less
% trivial
u   = tanh( (yu-pi/2)/d) .* (yu<=pi) + tanh( (3*pi/2 - yu)/d) .* (yu>pi);
v   = e*sin(xv);
p   = zeros(Npx,Npy);

    
end