% body force
% in Finite Volume setting, so integrated

% % only run if unsteady force or at first iteration
if (force_unsteady == 1 || exist('n','var')==0) 

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