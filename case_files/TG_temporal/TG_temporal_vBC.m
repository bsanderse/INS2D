function v = TG_temporal_vBC(x,y,t,options)
% boundary conditions for v for TG

Re = options.fluid.Re;

v = cos(pi*x)*sin(pi*y)*exp(-2*(pi^2)*t/Re);

end