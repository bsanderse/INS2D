function [phi_inhom,R_inhom,P] = get_phi_inhom(phi_bc,options)
% P only required for debugging

Mbc = size(phi_bc,2);

% Y_M = zeros(options.grid.Np,Mbc);
% for jj = 1:Mbc
%     yBC = phi_bc(:,jj);
%     Y_M(:,jj) = get_yM(options,yBC);
% end
% %% testing
% F_M = options.discretization.F_M;
% norm(Y_M-F_M*phi_bc)
%%
F_M = options.discretization.F_M;
Y_M = F_M*phi_bc;
%%

L = options.discretization.A;
Gx   = options.discretization.Gx;
Gy   = options.discretization.Gy;
G = [Gx;Gy];

Om = options.grid.Om;
Om_inv = options.grid.Om_inv;
%     tilde_phi_inhom = Om_inv.*(G*(L\Y_M));
tilde_phi_inhom = -Om_inv.*(G*(L\Y_M)); %pfusch

[phi_inhom,R_inhom,P] = Om_orthonormalize(tilde_phi_inhom,options);

%testing
% M_h = options.discretization.M;
% F_M = options.discretization.F_M;
% norm(M_h*phi_inhom*R_inhom-F_M*phi_bc(:,P))
% norm(phi_inhom*R_inhom-tilde_phi_inhom(:,P))


%% testing
% F_M = options.discretization.F_M;
% norm(yM-F_M*yBC)