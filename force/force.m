function [Fx, Fy] = force(t,options)
% body force
% in Finite Volume setting, so integrated

force_unsteady = options.case.force_unsteady;

% % only run if unsteady force or at first iteration of steady 
if (force_unsteady == 1 || t == options.time.t_start) 

    file_name = [options.case.project '_force'];
    
    if (exist(file_name,'file'))
        % create function handle with name bodyforce
        bodyforce = str2func(file_name);    
        [Fx, Fy]  = bodyforce(t,options);
    else
        % 
        Fx = zeros(options.grid.Nu,1);
        Fy = zeros(options.grid.Nv,1);
        
        if (force_unsteady == 1)
            error(['Body force file ' file_name ' not available']);
        end
    end
    
else
    Fx = zeros(options.grid.Nu,1);
    Fy = zeros(options.grid.Nv,1);
end

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
