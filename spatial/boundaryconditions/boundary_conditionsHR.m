%% get BC type

file_name = 'BCtype';

if (exist(file_name,'file'))
    
    BCtype     = str2func(file_name);
    options.BC = BCtype();
%     options.BC = BCtype(j);
    
else
    
    error(['BCtype file ' file_name ' not available']);
    
end

% check types
normal_BC = {'dir','per','pres','mvp-obc'};
tangential_BC = {'dir','per','sym','mvp-obc'};
if ( max(strcmp(options.BC.u.left,normal_BC))==0)
    error('wrong BC for u-left');
end
if ( max(strcmp(options.BC.u.right,normal_BC))==0)
    error('wrong BC for u-right');
end
if ( max(strcmp(options.BC.u.low,tangential_BC))==0)
    error('wrong BC for u-low');
end
if ( max(strcmp(options.BC.u.up,tangential_BC))==0)
    error('wrong BC for u-up');
end
if ( max(strcmp(options.BC.v.left,tangential_BC))==0)
    error('wrong BC for v-left');
end
if ( max(strcmp(options.BC.v.right,tangential_BC))==0)
    error('wrong BC for v-right');
end
if ( max(strcmp(options.BC.v.low,normal_BC))==0)
    error('wrong BC for v-low');
end
if ( max(strcmp(options.BC.v.up,normal_BC))==0)
    error('wrong BC for v-up');
end

%% check whether mvp-obc in u iff mvp-obc in v
if xor(strcmp(options.BC.u.left,'mvp-obc'),strcmp(options.BC.v.left,'mvp-obc'))
    error('BC mvp-obc must be set for both u and v')
end
if xor(strcmp(options.BC.u.right,'mvp-obc'),strcmp(options.BC.v.right,'mvp-obc'))
    error('BC mvp-obc must be set for both u and v')
end
if xor(strcmp(options.BC.u.low,'mvp-obc'),strcmp(options.BC.v.low,'mvp-obc'))
    error('BC mvp-obc must be set for both u and v')
end
if xor(strcmp(options.BC.u.up,'mvp-obc'),strcmp(options.BC.v.up,'mvp-obc'))
    error('BC mvp-obc must be set for both u and v')
end

%% set BC functions

% values set below can be either Dirichlet or Neumann value,
% depending on BC set above. in case of Neumann (symmetry, pressure)
% one uses normally zero gradient

% values should either be scalars or vectors
% ALL VALUES (u,v,p,k,e) are defined at x,y locations,
% i.e. the corners of pressure volumes, so they cover the entire domain
% including corners

% file_name = [options.case.project '_uBC'];
file_name = 'uBC';
if (exist(file_name,'file'))
    % create function handle with name uBC
    uBC = str2func(file_name);
else
    error(['BCtype file ' file_name ' not available']);
end

% uLo      = uBC(x,y(1),t,options);
% uUp      = uBC(x,y(end),t,options);
% uLe      = uBC(x(1),y,t,options);
% uRi      = uBC(x(end),y,t,options);


% file_name = [options.case.project '_vBC'];
file_name = 'vBC';
if (exist(file_name,'file'))
    % create function handle with name vBC
    vBC   = str2func(file_name);
else
    error(['BCtype file ' file_name ' not available']);
end
% vLo      = vBC(x,y(1),t,options);
% vUp      = vBC(x,y(end),t,options);
% vLe      = vBC(x(1),y,t,options);
% vRi      = vBC(x(end),y,t,options);

% time derivative of boundary conditions is used for high accuracy of
% pressure
% file_name = [options.case.project '_dudtBC'];
file_name = 'dudtBC';
if (exist(file_name,'file'))
    % create function handle with name dudtBC
    dudtBC = str2func(file_name);
else
    if (options.BC.BC_unsteady==1 && options.solversettings.p_add_solve==1)
        error(['Unsteady boundary conditions with additional Poisson solve requires time derivatives, but file ' file_name ' is not available']);
    end
end

% time derivative of boundary conditions is used for high accuracy of
% pressure
% file_name = [options.case.project '_dvdtBC'];
file_name = 'dvdtBC';
if (exist(file_name,'file'))
    % create function handle with name dudtBC
    dvdtBC = str2func(file_name);
else
    if (options.BC.BC_unsteady==1 && options.solversettings.p_add_solve==1)
        error(['Unsteady boundary conditions with additional Poisson solve requires time derivatives, but file ' file_name ' is not available']);
    end
end

%% pressure
% pressure BC is only used when at the corresponding boundary
% 'pres' is specified
% note: in the future, this should become an inputfile, e.g. casefile_pBC.m
p_inf    = 0;
pLe      = p_inf;
pRi      = p_inf;
pLo      = p_inf;
pUp      = p_inf;

options.BC.pLe = pLe;
options.BC.pRi = pRi;
options.BC.pLo = pLo;
options.BC.pUp = pUp;

%% k-eps values
% note: in the future, this should become an inputfile, e.g. casefile_kepsBC.m

switch visc
    
    case 'keps'
                
        kLo      = 0;
        kUp      = 0;
        kLe      = 0;
        kRi      = 0;
        
        eLo      = 0;
        eUp      = 0;
        eLe      = 0;
        eRi      = 0;
        
        options.BC.kLe = kLe;
        options.BC.kRi = kRi;
        options.BC.kLo = kLo;
        options.BC.kUp = kUp;
        
        options.BC.eLe = eLe;
        options.BC.eRi = eRi;
        options.BC.eLo = eLo;
        options.BC.eUp = eUp;
        
end


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
