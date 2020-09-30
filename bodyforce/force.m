function [Fx, Fy, dFx, dFy] = force(V,t,options,getJacobian)
% body force in momentum equations
% in Finite Volume setting, so integrated
% dFx, dFy are the Jacobians dFx/dV and dFy/dV

force_unsteady = options.case.force_unsteady;


if (options.force.isforce == 1)
    % create function handle with name bodyforce
    [Fx, Fy, dFx, dFy]  = options.force.bodyforce(V,t,options,getJacobian);
else
    % 
    Nu = options.grid.Nu;
    Nv = options.grid.Nv;
    Fx = zeros(Nu,1);
    Fy = zeros(Nv,1);
    dFx = spalloc(Nu,Nu+Nv,0);
    dFy = spalloc(Nv,Nu+Nv,0);    
    
    if (force_unsteady == 1)
        error(['Body force file ' file_name ' not available']);
    end
end
    
% else
%     Fx = zeros(options.grid.Nu,1);
%     Fy = zeros(options.grid.Nv,1);
% end

%% actuator methods:
% force_airfoil_surfaceCp;
% force_airfoil_camberCp;    
% force_airfoil_chordCp;
% force_airfoil_chordCL;

% force_actuatordisk;
% force_rotatedactuatordisk;
% force_rotatingactuatordisk;
% force_point;

% airfoil
% force_airfoil;

%% channel flow
% force_channel;


%% MMS
% force_MMS; 


%% IBM
% gcibm;


%% LDC - Shih
% force_LDC_steady;

%% Taylor vortex without pressure
% force_Taylor;

% force_turbines;

%% fictitious force
% Fx = Fx(:)+sin(t)*Omu;
% Fy = zeros(Nv,1);
% Fx = -sin(t)*Omu;
% Fy = zeros(Nv,1);

end
