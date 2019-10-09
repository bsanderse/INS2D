function BC = shear_layer_NOACK_ROM_BCtype

    BC.BC_unsteady  = 1;

    % 'dir' : inflow, wall
    % 'sym' : symmetry
    % 'pres': pressure (outflow)
    % 'per' : periodic

    % left/right: x-direction
    % low/up: y-direction

    BC.u.left  = 'dir';   % valid options: dir, per, pres 
    BC.u.right = 'pres';  % valid options: dir, per, pres
    BC.u.low   = 'sym';   % valid options: dir, per, sym
    BC.u.up    = 'sym';   % valid options: dir, per, sym

    BC.v.left  = 'dir';   % valid options: dir, per, sym
    BC.v.right = 'sym';   % valid options: dir, per, sym
    BC.v.low   = 'dir';   % valid options: dir, per, pres
    BC.v.up    = 'dir';   % valid options: dir, per, pres

