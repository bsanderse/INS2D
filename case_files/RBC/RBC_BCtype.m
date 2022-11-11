function BC = RBC_BCtype

    BC.BC_unsteady  = 0;

    % 'dir' : inflow, wall
    % 'sym' : symmetry
    % 'pres': pressure (outflow)
    % 'per' : periodic

    % left/right: x-direction
    % low/up: y-direction

    BC.u.left  = 'dir';   % valid options: dir, per, pres 
    BC.u.right = 'dir';  % valid options: dir, per, pres
    BC.u.low   = 'per';   % valid options: dir, per, sym
    BC.u.up    = 'per';   % valid options: dir, per, sym

    BC.v.left  = 'dir';   % valid options: dir, per, sym
    BC.v.right = 'dir';   % valid options: dir, per, sym
    BC.v.low   = 'per';   % valid options: dir, per, pres
    BC.v.up    = 'per';   % valid options: dir, per, pres

    BC.T.left  = 'sym';   % valid options: dir, sym, per
    BC.T.right = 'sym';   % valid options: dir, sym, per
    BC.T.low   = 'per';   % valid options: dir, sym, per
    BC.T.up    = 'per';   % valid options: dir, sym, per

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
