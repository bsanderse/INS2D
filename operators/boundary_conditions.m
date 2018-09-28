%% get BC type 

file_name = [options.case.project '_BCtype'];

if (exist(file_name,'file'))
    
    BCtype     = str2func(file_name);    
    options.BC = BCtype();

else
    
    error(['BCtype file ' file_name ' not available']);
    
end

%% get BC values

% values set below can be either Dirichlet or Neumann value, 
% depending on BC set above. in case of Neumann (symmetry, pressure) 
% one uses normally zero gradient

% values should either be scalars or vectors
% ALL VALUES (u,v,p,k,e) are defined at x,y locations, 
% i.e. the corners of pressure volumes, so they cover the entire domain
% including corners

file_name = [options.case.project '_uBC'];
if (exist(file_name,'file'))
    uBC  = str2func(file_name);    
else
    error(['BCtype file ' file_name ' not available']);
end

uLo      = uBC(x,y(1),t,options);
uUp      = uBC(x,y(end),t,options); 
uLe      = uBC(x(1),y,t,options);
uRi      = uBC(x(end),y,t,options);


file_name = [options.case.project '_vBC'];
if (exist(file_name,'file'))
    vBC   = str2func(file_name);    
else
    error(['BCtype file ' file_name ' not available']);
end
vLo      = vBC(x,y(1),t,options); 
vUp      = vBC(x,y(end),t,options);
vLe      = vBC(x(1),y,t,options); 
vRi      = vBC(x(end),y,t,options);


%% pressure
% pressure BC is only used when at the corresponding boundary 
% 'pres' is specified
p_inf    = 0;
pLe      = p_inf;
%     lambda   = Re/2-sqrt(Re^2/4+4*pi^2);   
pRi      = p_inf; %-0.5*exp(2*lambda*x2)+(lambda/Re)*exp(lambda*x2)*cos(2*pi*y);
pLo      = p_inf;
pUp      = p_inf;


%% k-eps values

kLo      = 0;
kUp      = 0; 
kLe      = 0;%(u_fr/U_ref)^2/sqrt(Cmu); 
kRi      = 0;

eLo      = 0;
eUp      = 0;
eLe      = 0;%kappa^2./(log((1+z0_nondim)/z0_nondim)^3*(y+z0_nondim)); 
eRi      = 0;
   
    
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
