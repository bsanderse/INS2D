%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% boundary conditions
    BC_unsteady  = 0;

    % 'dir' : inflow, wall
    % 'sym' : symmetry
    % 'pres': pressure (outflow)
    % 'per' : periodic

    % left/right: x-direction
    % low/up: y-direction

    BC.u.left  = 'dir';   % valid options: dir, per, pres 
    BC.u.right = 'dir';  % valid options: dir, per, pres
    BC.u.low   = 'dir';   % valid options: dir, per, sym
    BC.u.up    = 'dir';   % valid options: dir, per, sym

    BC.v.left  = 'dir';   % valid options: dir, per, sym
    BC.v.right = 'dir';   % valid options: dir, per, sym
    BC.v.low   = 'dir';   % valid options: dir, per, pres
    BC.v.up    = 'dir';   % valid options: dir, per, pres

    BC.k.left  = 'dir';   % valid options: dir, sym, per
    BC.k.right = 'dir';   % valid options: dir, sym, per
    BC.k.low   = 'dir';   % valid options: dir, sym, per
    BC.k.up    = 'dir';   % valid options: dir, sym, per

    BC.e.left  = 'dir';   % valid options: dir, sym, per
    BC.e.right = 'dir';   % valid options: dir, sym, per
    BC.e.low   = 'dir';   % valid options: dir, sym, per
    BC.e.up    = 'dir';   % valid options: dir, sym, per    
    
    % values set below can be either Dirichlet or Neumann value, 
    % depending on BC set above. in case of Neumann (symmetry, pressure) 
    % one uses normally zero gradient

    % values should either be scalars or vectors
    % ALL VALUES (u,v,p,k,e) are defined at x,y locations, 
    % i.e. the corners of pressure volumes, so they cover the entire domain
    % including corners
    uLo      = uBC(x,y(1),t,Re); %-sin(pi*x)*cos(pi*y(1))*exp(-2*pi^2*t/Re);
    uUp      = uBC(x,y(end),t,Re); %-sin(pi*x)*cos(pi*y(end))*exp(-2*pi^2*t/Re);%16*(x.^4-2*x.^3+x.^2)*cos(t); 
    uLe      = uBC(x(1),y,t,Re); %-sin(pi*x(1))*cos(pi*y)*exp(-2*pi^2*t/Re);
    uRi      = uBC(x(end),y,t,Re); %-sin(pi*x(end))*cos(pi*y)*exp(-2*pi^2*t/Re);

    vLo      = vBC(x,y(1),t,Re); %cos(pi*x)*sin(pi*y(1))*exp(-2*pi^2*t/Re);
    vUp      = vBC(x,y(end),t,Re); %cos(pi*x)*sin(pi*y(end))*exp(-2*pi^2*t/Re);
    vLe      = vBC(x(1),y,t,Re); %cos(pi*x(1))*sin(pi*y)*exp(-2*pi^2*t/Re);
    vRi      = vBC(x(end),y,t,Re); %cos(pi*x(end))*sin(pi*y)*exp(-2*pi^2*t/Re);
   
    % pressure BC is only used when at the corresponding boundary 
    % 'pres' is specified
    p_inf    = 0;
    pLe      = p_inf;
%     lambda   = Re/2-sqrt(Re^2/4+4*pi^2);   
    pRi      = p_inf; %-0.5*exp(2*lambda*x2)+(lambda/Re)*exp(lambda*x2)*cos(2*pi*y);
    pLo      = p_inf;
    pUp      = p_inf;

    kLo      = 0;
    kUp      = 0; 
    kLe      = 0;%(u_fr/U_ref)^2/sqrt(Cmu); 
    kRi      = 0;

    eLo      = 0;
    eUp      = 0;
    eLe      = 0;%kappa^2./(log((1+z0_nondim)/z0_nondim)^3*(y+z0_nondim)); 
    eRi      = 0;
   
    
    % Neumann BC used to extrapolate values to the boundary
    % change only in case of periodic to 'per', otherwise leave at 'sym'
%     BC.nu.left  = 'sym';   % 
%     BC.nu.right = 'sym';   % 
%     BC.nu.low   = 'sym';   % 
%     BC.nu.up    = 'sym';   % 
%     BC.nu.back  = 'sym';   % 
%     BC.nu.front = 'sym';   % 
%     nuLe        = 0;
%     nuRi        = 0;
%     nuLo        = 0;
%     nuUp        = 0;
%     nuBa        = 0;
%     nuFr        = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
