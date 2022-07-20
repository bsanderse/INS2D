function [phi_inhom,R_inhom] = get_phi_inhom(phi_bc,options)

Mbc = size(phi_bc,2);

for jj = 1:Mbc
    yBC = phi_bc(:,jj);
    Y_M(:,jj) = get_yM(options,yBC);
end
L = options.discretization.A;
Gx   = options.discretization.Gx;
Gy   = options.discretization.Gy;
G = [Gx;Gy];

Om = options.grid.Om;
Om_inv = options.grid.Om_inv;
%     tilde_phi_inhom = Om_inv.*(G*(L\Y_M));
tilde_phi_inhom = -Om_inv.*(G*(L\Y_M)); %pfusch

[Q_inhom,R_inhom] = qr(sqrt(Om).*tilde_phi_inhom,0); 
M_inhom = rank(tilde_phi_inhom);
Q_1_inhom = -Q_inhom(:,1:M_inhom);
R_inhom = -R_inhom(1:M_inhom,:);
phi_inhom = sqrt(Om_inv).*Q_1_inhom;

options.rom.phi_inhom = phi_inhom;
options.rom.R_inhom = R_inhom;