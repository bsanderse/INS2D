V_total = [uh_total vh_total]'; %[uh_total; vh_total];
snapshots_V_total = [test_snapshots.uh_total(snapshot_indx,:) test_snapshots.vh_total(snapshot_indx,:)]';
% error_V = V_total - snapshots_V_total;

M_h = options.discretization.M;
F_M = options.discretization.F_M;
phi_bc = options.rom.phi_bc;
Mbc = options.rom.Mbc;

% if j==1  % has to be done for each phi_bc separately
    no_ts = numel(t_vec);
    A_bc = zeros(Mbc,no_ts);
    for ii = 1:no_ts
        t_now = t_vec(ii);
        A_bc(:,ii) = get_a_bc(t_now,options);
    end
% end

% mass_viola = M_h*error_V;
y_Ms = F_M*phi_bc*A_bc;
% mass_viola = M_h*V_total - y_Ms;
mass_viola = M_h*V_total + y_Ms; % botch
% mass_viola_FOM = M_h*snapshots_V_total - y_Ms;
mass_viola_FOM = M_h*snapshots_V_total + y_Ms; % botch


Np = options.grid.Np;
one_vec = ones(Np,1);

mass_viola_L2 = weightedL2norm(mass_viola, one_vec);
% sum_mass_viola_L2(j) = sum(mass_viola_L2)
mass_viola_L2_FOM = weightedL2norm(mass_viola_FOM, one_vec);

disp('approx mass viola')
figure(112)
% plot(t_vec,mass_viola_L2,color,'displayname',"ROM M="+M+suffix);
semilogy(t_vec,mass_viola_L2,color,'displayname',"ROM M="+M+suffix);
sum2_mass_viola_L2(j) = sum(mass_viola_L2)
hold on
% semilogy(t_vec,mass_viola_L2_FOM,color2,'displayname',"FOM");

% if j==Nsim
%     semilogy(t_vec,snapshots.maxdiv(snapshot_indx),'displayname',"FOM");
% end

set(gca,'Yscale','log');
xlabel('t')
ylabel('approx mass equation violation')
legend('show','NumColumns',2,'Orientation','horizontal')

% figure(113)
% % plot(t_vec,dot(mass_viola_L2,error_p),color,'displayname',"pressure term
% % error ROM M="+M+suffix); %wrong
% % ROM_press_terms = dot(M_h*V_total,p_total'); % wrong: should take FOM mass eq RHS
% ROM_press_terms = dt*dot(M_h*snapshots_V_total,p_total'); % wrong: should take FOM mass eq RHS
% FOM_press_terms = dt*dot(M_h*snapshots_V_total,snapshots_p_total);
% 
% diff_pressure_terms = abs(ROM_press_terms-FOM_press_terms);
% plot(t_vec,diff_pressure_terms,color2,'displayname',"pressure term error ROM M="+M+suffix);
% 
% legend('show')
