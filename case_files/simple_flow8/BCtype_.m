function BC = BCtype_

    BC.BC_unsteady  = 0;

    % 'dir' : inflow, wall
    % 'sym' : symmetry
    % 'pres': pressure (outflow)
    % 'per' : periodic

    % left/right: x-direction
    % low/up: y-direction

    BC.u.left  = 'dir';   % valid options: dir, per, pres 
    BC.u.right = 'mvp-obc';  % valid options: dir, per, pres
%     BC.u.low   = 'mvp-obc';   % valid options: dir, per, sym
%     BC.u.up    = 'mvp-obc';   % valid options: dir, per, sym
    
%     BC.u.right = 'pres';  % valid options: dir, per, pres
    BC.u.low   = 'sym';   % valid options: dir, per, sym
    BC.u.up    = 'sym';   % valid options: dir, per, sym
    
%     BC.u.low   = 'dir';   % valid options: dir, per, sym
%     BC.u.up    = 'dir';   % valid options: dir, per, sym

    BC.v.left  = 'dir';   % valid options: dir, per, sym
    BC.v.right = 'mvp-obc';   % valid options: dir, per, sym
%     BC.v.low   = 'mvp-obc';   % valid options: dir, per, pres
%     BC.v.up    = 'mvp-obc';   % valid options: dir, per, pres

%     BC.v.right = 'sym';   % valid options: dir, per, sym
    BC.v.low   = 'pres';   % valid options: dir, per, pres
    BC.v.up    = 'pres';   % valid options: dir, per, pres

%     BC.v.low   = 'dir';   % valid options: dir, per, pres
%     BC.v.up    = 'dir';   % valid options: dir, per, pres

%     BC.k.left  = 'dir';   % valid options: dir, sym, per
%     BC.k.right = 'dir';   % valid options: dir, sym, per
%     BC.k.low   = 'dir';   % valid options: dir, sym, per
%     BC.k.up    = 'dir';   % valid options: dir, sym, per
% 
%     BC.e.left  = 'dir';   % valid options: dir, sym, per
%     BC.e.right = 'dir';   % valid options: dir, sym, per
%     BC.e.low   = 'dir';   % valid options: dir, sym, per
%     BC.e.up    = 'dir';   % valid options: dir, sym, per    
%     
    % values set below can be either Dirichlet or Neumann value, 
    % depending on BC set above. in case of Neumann (symmetry, pressure) 
    % one uses normally zero gradient


   
    
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

%     BC.gO = @(u) u.^2; % non-negative function, dissipation model %% maybe move to parameters
%     BC.gO = @(u) abs(u); % non-negative function, dissipation model %% maybe move to parameters
%     BC.gO = @(u) 1; % non-negative function, dissipation model %% maybe move to parameters
    BC.gO = @(u) 1; % non-negative function, dissipation model %% maybe move to parameters
%     BC.gO = @(u) 0; % non-negative function, dissipation model %% maybe move to parameters
% debugged Jacobian FOM data is with gO = 1!!!
    BC.dgO = @(u) 0; % non-negative function, dissipation model %% maybe move to parameters
%     BC.dgO = @(u) sign(u); % non-negative function, dissipation model %% maybe move to parameters
%     BC.dgO = @(u) 2*u; % non-negative function, dissipation model %% maybe move to parameters
    
    BC.gO_type = 1; % 0: gO=0, 1: gO=const, 2: gO more complex -> DEIM required
%     BC.gO_type = 0; % 0: gO=0, 1: gO=const, 2: gO more complex -> DEIM required

    
    BC.gO2string = func2str(BC.gO);
    BC.dgO2string = func2str(BC.dgO);
