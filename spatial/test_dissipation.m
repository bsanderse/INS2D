% test diffusion -> dissipation

% diffusion given by
% ((u_{i+1} - u_{i})/dx - (u_{i} - u_{i-1})/dx)/dx = (u_{i+1} - 2*u_{i}-
% u_{i-1})/dx^2

% we leave out the scaling with dx^2 below, it should be included in all
% terms
% note that in a FV setting these terms would be "multiplied" by dx*dy,
% leading effectively to a factor dy/dx in all terms

% upon product with u_{i} we get
% u_i*((u_{i+1} - u_{i}) - (u_{i} - u_{i-1})) = 
%  0.5*(u_{i+1} + u_{i})*(u_{i+1} - u_{i}) - 0.5*(u_{i+1} - u_{i})^2
% -0.5*(u_{i} + u_{i-1})*(u_{i} - u_{i-1}) - 0.5*(u_{i} - u_{i-1})^2
 


N = 6;
L = 1;


%% periodic case

e = ones(N,1);
x = linspace(0,L,N)';

u_per = sin(2*pi*x.^2/L);

Q_per = spdiags([-e e],[0 1],N,N); % take "right differences"
Q_per(end,1) = 1;
% note -Q_per' takes "left" differences
D_per = -Q_per'*Q_per; 
Ar_per = spdiags([0.5*e 0.5*e],[0 1],N,N); % take "right averages"
Ar_per(end,1) = 0.5;
Al_per = spdiags([0.5*e 0.5*e],[-1 0],N,N); % take "left" averages but without factor 0.5
Al_per(1,end) = 0.5;

% test u_i*((u_{i+1} - u_{i}) = 0.5*(u_{i+1} + u_{i})*(u_{i+1} - u_{i})
u_per.*(Q_per*u_per);
(Ar_per*u_per).*(Q_per*u_per) - 0.5*(Q_per*u_per).^2;

% now test full expression
u_per.*(D_per*u_per);
(-Q_per')*((Ar_per*u_per).*(Q_per*u_per)) - Al_per*((Q_per*u_per).^2);

%
sum(u_per.*(D_per*u_per))
- sum(Al_per*((Q_per*u_per).^2))

%% dirichlet case, points aligned with boundaries

e = ones(N+1,1);
x = linspace(0,L,N+1)'; 
% such that we have N-1 inner points (N 'pressure volumes')
x_dir = x(2:end-1);

% solution that satisfies homogeneous boundary conditions
% only defined at "inner" points
u_dir = x_dir.*(L-x_dir).^3;

% general i
% u_i*((u_{i+1} - u_{i}) - (u_{i} - u_{i-1})) = 
%  0.5*(u_{i+1} + u_{i})*(u_{i+1} - u_{i}) - 0.5*(u_{i+1} - u_{i})^2
% -0.5*(u_{i} + u_{i-1})*(u_{i} - u_{i-1}) - 0.5*(u_{i} - u_{i-1})^2

% i=1:
% u_1*((u_{2} - u_{1}) - (u_{1} - u_{b})) = 
%  0.5*(u_{2} + u_{1})*(u_{2} - u_{1}) - 0.5*(u_{2} - u_{1})^2
% -0.5*(u_{1} + u_{b})*(u_{1} - u_{b}) - 0.5*(u_{1} - u_{b})^2
% = (u_b=0)
%  0.5*(u_{2} + u_{1})*(u_{2} - u_{1}) - 0.5*(u_{2} - u_{1})^2
% -0.5*u_{1}^2 - 0.5*u_{1}^2
% (note that the term -0.5*(u_{1} + u_{b})*(u_{1} - u_{b}) does not
% evaluate to 0 for u_b=0)

% du/dx = Q_dir*u, with zero BC; this gives (u_1 - u_b), (u_2 - u_1), etc.
Q_dir = spdiags([-e e],[-1 0],N,N-1); 
D_dir = -Q_dir'*Q_dir;

% u_avg = Ar_dir *u, with zero BC; this gives (u_1 + u_b)/2, (u_2 + u_1)/2,
% etc.
Ar_dir = spdiags([0.5*e 0.5*e],[-1 0],N,N-1); 

Al_dir = spdiags([0.5*e 0.5*e],[0 1],N-1,N); 

% local expressions
u_dir.*(D_dir*u_dir);
(-Q_dir')*((Ar_dir*u_dir).*(Q_dir*u_dir)) - Al_dir*((Q_dir*u_dir).^2);

% test whether the conservative terms cancel upon summation
sum(u_dir.*(D_dir*u_dir))
- sum(Al_dir*(Q_dir*u_dir).^2)

% this is not the case because -0.5*(u_{1} + u_{b})*(u_{1} - u_{b}) does not
% evaluate to 0 for u_b=0, instead it evaluates to -0.5*u_{1}^2
% we therefore adapt the dissipation term to read:
- sum(Al_dir*(Q_dir*u_dir).^2) - 0.5*u_dir(1)^2 - 0.5*u_dir(end)^2
% and now they're equal to u'*(D*u)


%% dirichlet case, points not aligned with boundaries

e = ones(N+1,1);
x = linspace(0,L,N+2)'; 
% such that we have N inner points
x_dir = x(2:end-1);

% solution that satisfies homogeneous boundary conditions
% only defined at "inner" points
u_dir = x_dir.*(L-x_dir).^3;

% general i
% u_i*((u_{i+1} - u_{i}) - (u_{i} - u_{i-1})) = 
%  0.5*(u_{i+1} + u_{i})*(u_{i+1} - u_{i}) - 0.5*(u_{i+1} - u_{i})^2
% -0.5*(u_{i} + u_{i-1})*(u_{i} - u_{i-1}) - 0.5*(u_{i} - u_{i-1})^2

% i=1: u_0 + u_1 = 2*u_b, so u_0 = 2*u_b - u_1, and u_1 - u_0 = 2*u_1 - 2*ub
% u_1*((u_{2} - u_{1}) - (u_{1} - u_{0})) = 
%  0.5*(u_{2} + u_{1})*(u_{2} - u_{1}) - 0.5*(u_{2} - u_{1})^2
% -0.5*(u_{1} + u_{0})*(u_{1} - u_{0}) - 0.5*(u_{1} - u_{0})^2 =
%  0.5*(u_{2} + u_{1})*(u_{2} - u_{1}) - 0.5*(u_{2} - u_{1})^2
%               -2*u_{b}*(u_{1} - u_{b}) - 0.5*(2*u_{1} - 2*u_{b})^2
% = (u_b=0)
%  0.5*(u_{2} + u_{1})*(u_{2} - u_{1}) - 0.5*(u_{2} - u_{1})^2
% -2*(u_{1} - u_{b})^2
% (note that the term (u_{1} + u_{0})*(u_{1} - u_{0}) DOES evaluate to 0
% for u_b=0, in contrast to the aligned case

% du/dx = Q_dir*u, with zero BC; this gives (u_1 - u_0)=2*(u_1-u_b), (u_2 - u_1), etc.
Q_dir = spdiags([-e e],[-1 0],N+1,N); 
Q_dir(1,1) = 2;
Q_dir(end,end) = -2;
% take derivative back to u points
Q_dir2 = spdiags([-e e],[0 1],N,N+1); 
D_dir = Q_dir2*Q_dir;

% u_avg = Ar_dir *u, with zero BC; this gives (u_1 + u_0)/2 = u_b, (u_2 + u_1)/2,
% etc.
Ar_dir = spdiags([0.5*e 0.5*e],[-1 0],N+1,N); 
Ar_dir(1,1) = 0;
Ar_dir(end,end) = 0;

Al_dir = spdiags([0.5*e 0.5*e],[0 1],N,N+1); 

% local expressions
u_dir.*(D_dir*u_dir)
Q_dir2*((Ar_dir*u_dir).*(Q_dir*u_dir)) - Al_dir*((Q_dir*u_dir).^2)

% test whether the conservative terms cancel upon summation by comparing
% full expression to dissipation only
sum(u_dir.*(D_dir*u_dir))
- sum(Al_dir*(Q_dir*u_dir).^2)

% this is indeed the case! hurray

