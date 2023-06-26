V_total = [uh_total vh_total]'; %[uh_total; vh_total];
snapshots_V_total = [test_snapshots.uh_total(snapshot_indx,:) test_snapshots.vh_total(snapshot_indx,:)]';
error_V = V_total - snapshots_V_total;

M_h = options.discretization.M;

mass_viola = M_h*error_V;

Np = options.grid.Np;
one_vec = ones(Np,1);

mass_viola_L2 = weightedL2norm(mass_viola, one_vec);

disp('mass viola')
figure(111)
% plot(t_vec,mass_viola_L2,color,'displayname',"ROM M="+M+suffix);
% semilogy(t_vec,mass_viola_L2,color,'displayname',"ROM M="+M+suffix);
% semilogy(t_vec,mass_viola_L2,'color',color,'displayname',name);
semilogy(t_vec,mass_viola_L2,'color',color,'linestyle',linestyle,'linewidth', linewidth,'displayname',name);
sum_mass_viola_L2(j) = sum(mass_viola_L2)
hold on

if j==Nsim
    semilogy(t_vec,snapshots.maxdiv(snapshot_indx),'displayname',"FOM", 'color','k', 'linewidth', linewidth);
    set(gca,'Yscale','log');
    xlabel('M')
    ylabel('orthonormality error')
    legend('show')
end

set(gca,'Yscale','log');
xlabel('t')
ylabel('mass equation violation')
% legend('show','NumColumns',2,'Orientation','horizontal')
% legend('show','NumColumns',3,'Orientation','horizontal')
legend('show','NumColumns',3,'Orientation','vertical')
    set(gcf, 'Position', [100, 100, 400, 300])
    grid on

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


