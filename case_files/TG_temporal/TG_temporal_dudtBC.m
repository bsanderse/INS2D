function dudt = TG_temporal_dudtBC(x,y,t,options)

Re = options.fluid.Re;


% Taylor:
dudt = (2*(pi^2)/Re)*sin(pi*x)*cos(pi*y)*exp(-2*(pi^2)*t/Re);


end