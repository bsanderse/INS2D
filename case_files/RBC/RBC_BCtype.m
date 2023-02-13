function BC = RBC_BCtype

    BC.BC_unsteady  = 0;

    % 'dir' : inflow, wall
    % 'sym' : symmetry
    % 'pres': pressure (outflow)
    % 'per' : periodic

    % left/right: x-direction
    % low/up: y-direction

    temp=0; %temp =1 fpr periodic case and 0 for adiabatic and noslip case
if temp==0 
    BC.u.left  = 'dir';   % valid options: dir, per, pres 
    BC.u.right = 'dir';  % valid options: dir, per, pres
    BC.u.low   = 'dir';   % valid options: dir, per, sym
    BC.u.up    = 'dir';   % valid options: dir, per, sym

    BC.v.left  = 'dir';   % valid options: dir, per, sym
    BC.v.right = 'dir';   % valid options: dir, per, sym
    BC.v.low   = 'dir';   % valid options: dir, per, pres
    BC.v.up    = 'dir';   % valid options: dir, per, pres

    BC.T.left  = 'sym';   % valid options: dir, sym, per
    BC.T.right = 'sym';   % valid options: dir, sym, per
    BC.T.low   = 'dir';   % valid options: dir, sym, per
    BC.T.up    = 'dir';   % valid options: dir, sym, per
elseif temp==1 
    BC.u.left  = 'per';   % valid options: dir, per, pres 
    BC.u.right = 'per';  % valid options: dir, per, pres
    BC.u.low   = 'dir';   % valid options: dir, per, sym
    BC.u.up    = 'dir';   % valid options: dir, per, sym

    BC.v.left  = 'per';   % valid options: dir, per, sym
    BC.v.right = 'per';   % valid options: dir, per, sym
    BC.v.low   = 'dir';   % valid options: dir, per, pres
    BC.v.up    = 'dir';   % valid options: dir, per, pres

    BC.T.left  = 'per';   % valid options: dir, sym, per
    BC.T.right = 'per';   % valid options: dir, sym, per
    BC.T.low   = 'dir';   % valid options: dir, sym, per
    BC.T.up    = 'dir';   % valid options: dir, sym, per
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
