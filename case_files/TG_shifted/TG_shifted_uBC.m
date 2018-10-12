function u = TG_shifted_uBC(x,y,t,options)
% boundary conditions for u for TG

Re = options.fluid.Re;

% Taylor:
u = -sin(pi*x).*cos(pi*y)*exp(-2*(pi^2)*t/Re);


end