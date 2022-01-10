function BC = shear_layer_ROM_BCtype

    BC.BC_unsteady  = 0;

    % 'dir' : inflow, wall
    % 'sym' : symmetry
    % 'pres': pressure (outflow)
    % 'per' : periodic

    % left/right: x-direction
    % low/up: y-direction

    BC.u.left  = 'per';%'mvp-obc';   % valid options: dir, per, pres 
    BC.u.right = 'per';  % valid options: dir, per, pres
    BC.u.low   = 'mvp-obc';   % valid options: dir, per, sym
    BC.u.up    = 'mvp-obc';   % valid options: dir, per, sym

    BC.v.left  = 'per';   % valid options: dir, per, sym
    BC.v.right = 'per';   % valid options: dir, per, sym
    BC.v.low   = 'mvp-obc';   % valid options: dir, per, pres
    BC.v.up    = 'mvp-obc';   % valid options: dir, per, pres

%     BC.gO = @(u) 1; % non-negative function, dissipation model %% maybe move to parameters
    BC.gO = @(u) 0; % non-negative function, dissipation model %% maybe move to parameters
