% approx rom tests

Bphi = B'*(Om.*phi_inhom);
BB = B'*(Om.*B);
phiphi = phi_inhom'*(Om.*phi_inhom);

actuator_unsteady_ROM_init_unsteady_Vbc;

ratio = phi_inhom./Vbc__;
plot(ratio) %problem visible! -> actually no problem

factor = ratio(1);
diff = phi_inhom-factor*Vbc__;
figure
plot(diff); % no problem

ratio2 = tilde_phi_inhom(:,1)./Vbc__;
plot(ratio2) %problem visible! -> actually no problem

factor2 = ratio2(1);
diff2 = tilde_phi_inhom(:,1)-factor2*Vbc__;
figure
plot(diff2); % no problem

norm(-G'*tilde_phi_inhom(:,1)-Y_M(:,1)); %no problem
norm(-G'*tilde_phi_inhom-Y_M) %no problem

norm(-G'*Vbc__-Y_M(:,1))

% Y_q = L\Y_M;
% norm(L*Y_q-Y_M)
% norm(G*L*Y_q-G*Y_M)
% norm(Om_inv.*(G*L*Y_q)-Om_inv.*(G*Y_M))
% norm(Om_inv.*(G*L*Y_q)-Om_inv.*(G*Y_M))
% norm(tilde_phi_inhom-Om_inv.*(G*Y_q))
% norm(Y_M+G'*(Om_inv.*G*(L\Y_M)))
% norm(Y_q-L\Y_M)
% norm(G*Y_q-G*(L\Y_M))
% norm(Om_inv.*(G*Y_q)-Om_inv.*(G*(L\Y_M)))


for jj = 1:numel(t_js)
    t_j = t_js(jj);
    options = set_bc_vectors(t_j,options);
    Y_M2(:,jj) = options.discretization.yM;
    nnz_YM2(:,jj) = nonzeros(options.discretization.yM);
end
(Y_M2~=0)==(Y_M(:,1)~=0);

nnz_YM = nonzeros(Y_M(:,1));

nnz_ratio = nnz_YM2./nnz_YM;
surf(nnz_ratio)
mean = sum(nnz_ratio,1)/size(nnz_ratio,1);
devia = nnz_ratio - mean;
max(max(abs(devia)))


% rratio = Y_M2./Y_M(:,1);
% mean = sum(rratio,1)/size(rratio,1);
% devia = rratio - mean;


%%
actuator_unsteady_ROM_init_unsteady_Vbc;


t= 0;
Vbc2 = -options.rom.Vbc0*options.rom.abc(t);
Vbc3 = phi_inhom*get_a_inhom(t,options);
Vbc4 = tilde_phi_inhom*get_a_bc(t,options);

norm(Vbc2-Vbc3)
figure
plot(Vbc2-Vbc3)
% hold on
figure
plot(Vbc2-Vbc4)

psi = get_streamfunction(Vbc2,t,options);
% [up,vp,qp] = get_velocity(Vbc2,t,options);
figure(1)
% set(gcf,'color','w');
Axes(1) = axes;
contour(Axes(1),x(2:end-1),y(2:end-1),reshape(psi,Nx-1,Ny-1)','k'); %labels,'LineWidth',1);

% V_ = Vbc2;
V_ = Vbc3;
% V_ = Vbc4;
V_ = phi_inhom/norm(sqrt(Om).*phi_inhom);
V_ = Vbc__/norm(sqrt(Om).*Vbc__);
psi = get_streamfunction(V_,t,options);

[up,vp,qp] = get_velocity(V_,t,options);
% list = 20;
figure
set(gcf,'color','w');

Axes(1) = axes;
contour(Axes(1),x(2:end-1),y(2:end-1),reshape(psi,Nx-1,Ny-1)','k'); %labels,'LineWidth',1);
% colorbar
axis equal
axis([x1 x2 y1 y2]);
hold off

Axes(2) = axes;
% pcolor(xp,yp,qp')
% list = linspace(0.6,1.1,20);
list = linspace(-1,1,50);

[~,c]=contour(Axes(2),xp,yp,qp',list);
c.LineWidth = 1;
axis equal
axis([x1 x2 y1 y2]);
colorbar('Location','east')

%%
Vbc_2 = -Vbc__/norm(sqrt(Om).*Vbc__);
Vbc_3 = phi_inhom/norm(sqrt(Om).*phi_inhom);

norm(Vbc_2-Vbc_3)
