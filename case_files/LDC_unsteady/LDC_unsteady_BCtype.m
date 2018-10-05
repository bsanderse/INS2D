function BC = LDC_BCtype

    BC.BC_unsteady  = 0;

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
