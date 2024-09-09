function [u,v,p,options] = shear_layer_ROM_IC(t,options)
% initial velocity field Taylor-Green

yu = options.grid.yu;
xv = options.grid.xv;

Npx = options.grid.Npx;
Npy = options.grid.Npy;

%% original stuff
% note: we add 1 to u in order to make global momentum conservation less
% trivial
offset = 1;
d   = pi/15;
e   = 0.05;

u   = offset + tanh( (yu-pi/2)/d) .* (yu<=pi) + tanh( (3*pi/2 - yu)/d) .* (yu>pi);
v   = e*sin(xv);
p   = zeros(Npx,Npy);

%% opinf stuff

% % offset = -1;
% % d   = pi/15;
% % e   = 0.05;
% 
% rotate = options.case.IC_params(1);
% d   = options.case.IC_params(2);
% e   = options.case.IC_params(3);
% 
% % note: we add 1 to u in order to make global momentum conservation less
% % trivial
% u   = tanh( (yu-pi/2)/d) .* (yu<=pi) + tanh( (3*pi/2 - yu)/d) .* (yu>pi);
% v   = e*sin(xv);
% p   = zeros(Npx,Npy);
% 
% if rotate
%     % temp = u;
%     % u = v;
%     % v = temp;
% 
%     v   = tanh( (xv-pi/2)/d) .* (xv<=pi) + tanh( (3*pi/2 - xv)/d) .* (xv>pi);
%     u   = e*sin(yu);
% end
% 
% 
% % yv = options.grid.yu;
% % xu = options.grid.xv;
% 
% % v   = offset + tanh( (xv-pi/2)/d) .* (xv<=pi) + tanh( (3*pi/2 - xv)/d) .* (xv>pi);
% % u   = e*sin(yu);
% 
% disp("Divergence error of initial condition = " + num2str(norm(options.discretization.M*[u(:);v(:)])))
% 
% 
% 
% % botch
% % u = 0*u + 1;
% % v = 0*v + 1;
% 
% % % other botch: horizontal periodic
% % u = yu;
% % v = e*sin(xv);
% % 
% % % other botch: vertical periodic
% % u = e*sin(yu);
% % v = xv;
    
end