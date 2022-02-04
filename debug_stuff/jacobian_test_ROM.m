% jacobian test ROM
% [maxres,Fres,dFres] = F_ROM(R,~,t,options,getJacobian)

%%
    Conv_diag = options.grid.C;
    B = options.rom.B;

    M   = options.rom.M;
    E   = speye(M);

func3 = @(R) options.rom.obc_hom*kron(R,R);

Jac3 = options.rom.obc_hom*kron(E,R) + options.rom.obc_hom*kron(R,E);

% func = @(R) options.rom.obc_hom*kron(R,R);
% 
% Jac = options.rom.obc_hom*kron(E,R) + options.rom.obc_hom*kron(R,E);

%%
    M   = options.rom.M;
    E   = speye(M);
R_inhom = get_a_inhom(t,options);
R_bc    = get_a_bc(t,options);

Diff  = options.rom.Diff2;
DiffBC  = options.rom.DiffBC2;
yDiff = options.rom.yDiff2;

func0 = @(R) options.rom.C_hom2*kron(R,R) ...
     + options.rom.C_hom_inhom2*kron(R,R_inhom) ...
     + options.rom.C_hom_bc2*kron(R,R_bc) ...
     + options.rom.C_inhom2*kron(R_inhom,R_inhom) ...
     + options.rom.C_inhom_bc2*kron(R_inhom,R_bc) ...
     + options.rom.C_bc2*kron(R_bc,R_bc);
 
func = @(R) -func0(R)  +Diff*R + DiffBC*R_inhom + yDiff*R_bc;

% F00 = -func(R)+Diff*R + DiffBC*R_inhom + yDiff*R_bc;

    Jac = options.rom.C_hom2*(kron(E,R)+kron(R,E)) ...
        + options.rom.C_hom_inhom2*kron(E,R_inhom) ...
        + options.rom.C_hom_bc2*kron(E,R_bc);
% Jac00 = -Jac + Diff;
Jac = -Jac + Diff;

func = @(R) func(R) + func3(R);
Jac = Jac + Jac3;

%%
% Ct = options.force.Ct;
% options.force.Ct = 0;
% Re = options.fluid.Re;
% options.fluid.Re = 10^100;

%%
    [~,F0,Jac] = F_ROM(R,0,t,options,1);

%         F0 = func(R);
%         
%     [~,F00,Jac00] = F_ROM(R,0,t,options,1);
% norm(Jac00-Jac)
% norm(F0-F00)
%%
%     
%     VR = B*R;
%     F00 = B'*diag(VR)*(Conv_diag*VR);
%     norm(F0-F00) % = machine precision
% 
%     Jac00 = B'*(diag(Conv_diag*VR)+diag(VR)*Conv_diag);
%     norm(Jac*B'-Jac00) % = machine precision
    
%%

epsilons = 10.^(-1*(-4:8));
% epsilons = 10.^(-3);

errors = 0*epsilons;
for jj = 1:numel(epsilons)
epsilon = epsilons(jj);
Jac2 = zeros(numel(R));

for i = 1:numel(R)
    Ri = R;
    Ri(i) = Ri(i)+epsilon;
    
    [~,Fi,~] = F_ROM(Ri,0,t,options,1);

%     [~,Fii,~] = F_ROM(Ri,0,t,options,1);
%     Fi = func(Ri);
%     norm(Fii-Fi)
    
    dFi = (Fi-F0)/epsilon;
    Jac2(i,:) = dFi;
end

errorj = norm(Jac-Jac2');
errors(jj) = errorj;

end
figure(177)
loglog(epsilons,errors,'x')

1
% options.force.Ct = Ct;
% options.fluid.Re = Re;