function BC = shear_layer_ROM_BCtype
% function BC = shear_layer_ROM_BCtype(j)

    BC.BC_unsteady  = 0;

    % 'dir' : inflow, wall
    % 'sym' : symmetry
    % 'pres': pressure (outflow)
    % 'per' : periodic

    % left/right: x-direction
    % low/up: y-direction
% 
if j == 1
    BC.u.left  = 'per';%'mvp-obc';   % valid options: dir, per, pres 
    BC.u.right = 'per';  % valid options: dir, per, pres
    BC.u.low   = 'mvp-obc';   % valid options: dir, per, sym
    BC.u.up    = 'mvp-obc';   % valid options: dir, per, sym

    BC.v.left  = 'per';   % valid options: dir, per, sym
    BC.v.right = 'per';   % valid options: dir, per, sym
    BC.v.low   = 'mvp-obc';   % valid options: dir, per, pres
    BC.v.up    = 'mvp-obc';   % valid options: dir, per, pres
%%  
else
    BC.u.left  = 'pres';%'mvp-obc';   % valid options: dir, per, pres 
    BC.u.right = 'pres';  % valid options: dir, per, pres
    BC.u.low   = 'per';   % valid options: dir, per, sym
    BC.u.up    = 'per';   % valid options: dir, per, sym

    BC.v.left  = 'sym';   % valid options: dir, per, sym
    BC.v.right = 'sym';   % valid options: dir, per, sym
    BC.v.low   = 'per';   % valid options: dir, per, pres
    BC.v.up    = 'per';   % valid options: dir, per, pres
end
%%
%     BC.u.left  = 'per';%'mvp-obc';   % valid options: dir, per, pres 
%     BC.u.right = 'per';  % valid options: dir, per, pres
%     BC.u.low   = 'sym';   % valid options: dir, per, sym
%     BC.u.up    = 'sym';   % valid options: dir, per, sym
% 
%     BC.v.left  = 'per';   % valid options: dir, per, sym
%     BC.v.right = 'per';   % valid options: dir, per, sym
%     BC.v.low   = 'pres';   % valid options: dir, per, pres
%     BC.v.up    = 'pres';   % valid options: dir, per, pres
%%
    
%     BC.gO = @(u) 1; % non-negative function, dissipation model %% maybe move to parameters
    BC.gO = @(u) 0; % non-negative function, dissipation model %% maybe move to parameters

    BC.dgO = @(u) 0; % non-negative function, dissipation model %% maybe move to parameters
    
    BC.gO_type = 0; % 0: gO=0, 1: gO=const, 2: gO more complex -> DEIM required
