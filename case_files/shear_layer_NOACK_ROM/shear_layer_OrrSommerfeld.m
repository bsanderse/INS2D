function [ u,v,omega ] = shear_layer_OrrSommerfeld( alpha, t ,options )
%% Orr-Sommerfeld analysis for single-phase flow (2D)

%
% see http://www.chebfun.org/examples/ode-eig/OrrSommerfeld.html,
% http://www.chebfun.org/docs/guide/guide07.html, 
% and book of Schmid, Henningson, equation (3.14) or (3.19) and figure (3.14)
% 
addpath('/Users/sanderse/Dropbox/work/Programming/libs/chebfun-master/');
% addpath('/Users/sanderse/Dropbox/work/Programming/libs/');
% addpath('/Users/sanderse/Dropbox/work/Programming/libs/altmany-export_fig-cafc7c5/');

% close all
% clearvars
% clc

%% the analysis is based on the incompressible Navier-Stokes equations,
% which are linearized around a base state - tanh flow, see 
% Noack (JFM 2005): the need for a pressure-term representation of empirical
% Galerkin models 

% a fourth order equation for the perturbations v' can be derived, which
% reduces to an eigenvalue problem when assuming  the form 
% v' = v(y) exp(i(alpha*x-omega*t))
I     = sqrt(-1);


%% parameters
Re  = options.fluid.Re;

dom = [options.grid.y1,options.grid.y2];

% base flow: 
% U(y) = 0.5*(U1+U2) + 0.5*(U1-U2)*tanh(y/d), 
% U''(y) = 0.5*(U1-U2)*(2*tanh(y/d)*(tanh(y/d)^2 - 1))/(d^2)
U1 = options.fluid.U1;
U2 = options.fluid.U2;
d = options.fluid.d_layer;


choose_mode = 1; % with sorting based on largest imaginary part

neigs = 10;

%% set up chebfun operator (operator on functions)

A     = chebop(dom); % empty chebop on [-1,1]
A.op  = @(y,v) (  I*alpha*(0.5*(U1+U2) + 0.5*(U1-U2).*tanh(y/d)).*(diff(v,2)-alpha^2*v) + ...
                 -I*alpha*(0.5*(U1-U2)*(2*tanh(y/d).*(tanh(y/d).^2 - 1))/(d^2)).*v + ...
                 -(1/Re)*(diff(v,4)-2*(alpha^2)*diff(v,2)+(alpha^4)*v) );
             
B     = chebop(dom);
B.op  = @(y,v) (I*alpha)*(diff(v,2) - (alpha^2)*v); 

% homogeneous boundary conditions for v' and u' -> dv'/dy (only needed for A)
% at both sides:

if (strcmp(options.BC.v.low,'dir') && strcmp(options.BC.u.low,'dir')) 
        A.lbc = @(v) [v; diff(v,1)];     
elseif (strcmp(options.BC.v.low,'dir') && strcmp(options.BC.u.low,'sym'))
        A.lbc = @(v) [v; diff(v,2)];     
else
    error('this BC type not implemented in Orr-Sommerfeld');
end

if (strcmp(options.BC.v.up,'dir') && strcmp(options.BC.u.up,'dir')) 
        A.rbc = @(v) [v; diff(v,1)];     
elseif (strcmp(options.BC.v.up,'dir') && strcmp(options.BC.u.up,'sym'))
        A.rbc = @(v) [v; diff(v,2)];     
else
    error('this BC type not implemented in Orr-Sommerfeld');
end
% A.lbc = [0; 0];
% A.rbc = [0; 0];
% % or: 
% A.lbc = @(v) [v; diff(v,1)];     
% A.rbc = @(v) [v; diff(v,1)];     

% homogeneous boundary conditions for v' and du'/dy -> d2v'/dy2 (only needed for A)
% A.lbc = @(v) [v; diff(v,2)];     
% A.rbc = @(v) [v; diff(v,2)];     


% solve eigenvalue problem for c (phase speed) = omega/alpha
% A v = c B v, where the c is the same as in Schmid, eq (3.19)

% [V,D]  = eigs(A,B,60,'LR'); % LR: largest real part
[V,D]  = eigs(A,B,neigs,'LI'); % LI: largest imaginary part

% sort eigenvalues according to imaginary part
c        = diag(D);
[~,indx] = sort(imag(c),'desc');
% most unstable / least stable eigenmode:
c(indx(1));
c(indx(2));

maxe = max(real(c));

if (maxe > 0)
    disp('unstable');
else
    disp('stable');
end

% omega = c*alpha

omega = c(indx(choose_mode))*alpha;
disp(['omega=' num2str(omega)]);

% get perturbations as chebfuns

% v-component, perturbation
% convert chebmatrix to chebfun
v  = chebfun(V(:,indx(choose_mode)));

% u-component, perturbation
% u' = -1/(i*alpha) * dv/dy
% the chebfun can simply be differentiated
u  = (-1/(I*alpha)) * diff(v);



% contourf(real(reshape(u,options.grid.Nux_in,options.grid.Nuy_in)'))
% contourf(real(reshape(v,options.grid.Nvx_in,options.grid.Nvy_in)'))

end

