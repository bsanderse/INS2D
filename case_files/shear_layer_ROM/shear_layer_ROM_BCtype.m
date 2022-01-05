function BC = shear_layer_ROM_BCtype

    BC.BC_unsteady  = 0;

    % 'dir' : inflow, wall
    % 'sym' : symmetry
    % 'pres': pressure (outflow)
    % 'per' : periodic

    % left/right: x-direction
    % low/up: y-direction

    BC.u.left  = 'mvp-obc';%'mvp-obc';   % valid options: dir, per, pres 
    BC.u.right = 'mvp-obc';  % valid options: dir, per, pres
    BC.u.low   = 'dir';   % valid options: dir, per, sym
    BC.u.up    = 'dir';   % valid options: dir, per, sym

    BC.v.left  = 'mvp-obc';   % valid options: dir, per, sym
    BC.v.right = 'mvp-obc';   % valid options: dir, per, sym
    BC.v.low   = 'dir';   % valid options: dir, per, pres
    BC.v.up    = 'dir';   % valid options: dir, per, pres

%     BC.gO = @(u) 1; % non-negative function, dissipation model %% maybe move to parameters
    BC.gO = @(u) 0; % non-negative function, dissipation model %% maybe move to parameters
