
%%
% Ttest = rand(Npy,1);
% T = kron(Ttest,ones(Npx,1));
% T2D = rand(Npx,Npy);
% T = T2D(:);

Ttest = rand(Npx,1);
T = kron(ones(Npy,1),Ttest);

dTdy = options.discretization.STy*T+options.discretization.ySTy;
dTdyL = dTdy(1:Npx);
dTdyU = dTdy(end-Npx+1:end);

NusseltL = sum(-dTdyL.*hx) % integrate over lower plate
NusseltU = sum(-dTdyU.*hx) % integrate over upper plate


diffT    = diffusion_temperature(T,t,options,0);

epsilonT = thermal_dissipation(T,t,options);


% we want to show that
% int (T*d2T/dy2) dx dy = int T*(dT/dy) dx - int (dT/dy)^2 dx dy
% since T=0 on upper plate and T=1 on lower plate:
% T'*(DiffT + yT) = alfa4*NusseltL - epsilonT
T'*diffT/options.temp.alfa4
NusseltL - epsilonT/options.temp.alfa4
% -options.temp.alfa4*NusseltU - epsilonT



%% 
N=Npy;
dy = 1/N;
dx = 1/N;
y = (0.5:N-0.5)'/N;
e=ones(N,1);
DiffTest = spdiags([e -2*e e],[-1 0 1],N,N) * (dx/dy);
DiffTest(1,1) = -3*(dx/dy);
DiffTest(end,end) = -3*(dx/dy);
TLo = 1;
TUp = 0;
yTest= zeros(N,1);
yTest(1) = 2*TLo*(dx/dy);
yTest(N) = 2*TUp*(dx/dy);

% Ttest = 1-y.^5;
Ttest = rand(Npy,1);
% Ttest = T2D(:,1);

Ttest'*(DiffTest*Ttest + yTest)

dTdy = zeros(N+1,1);
dTdy(2:N)  = (Ttest(2:end) - Ttest(1:end-1))/dy;
dTdy(1)    = (Ttest(1) - TLo)/(0.5*dy);
dTdy(end)  = (TUp- Ttest(end))/(0.5*dy);
dTdy2      = dTdy.^2;
dTdy2(1)   = 0.5*dTdy2(1);
dTdy2(end) = 0.5*dTdy2(end);
eps = sum(dTdy2)*dx*dy;

-eps - TLo*(Ttest(1)-TLo)*dx/(0.5*dy)
% + TUp*(TUp-Ttest(end))*dx/(0.5*dy)

keyboard