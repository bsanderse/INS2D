V_total = [uh_total vh_total]'; %[uh_total; vh_total];
snapshots_V_total = [test_snapshots.uh_total(snapshot_indx,:) test_snapshots.vh_total(snapshot_indx,:)]';
error_V = V_total - snapshots_V_total;
% inf-norm
error_V_inf = max(abs(error_V),[],1);

% V_ref    = ones(size(snapshots_V_total));
V_ref    = snapshots_V_total;
figure(99)

V_2_ref  = weightedL2norm(V_ref,options.grid.Om);
V_2_ref_avg = sum(V_2_ref)/numel(V_2_ref)
error_V_2 = weightedL2norm(error_V,options.grid.Om)/V_2_ref_avg; ylabel('relative velocity error')

% error_V_2 = weightedL2norm(error_V,options.grid.Om); ylabel('velocity error')


%% wrong
% best possible approximation given the projection:
% note that the norm should be consistent with the
% optimization problem used in the SVD
% B = options.rom.B;
% Om = options.grid.Om;
% if options.verbosity.equivalence_cheat
%     V_best = B*(B'*(Om.*(snapshots_V_total)));
% else
%     V_best = B*(B'*(Om.*(snapshots_V_total)))+snapshots.Vbc;
% end
%% correct
B = options.rom.B;
Om = options.grid.Om;
% if options.verbosity.equivalence_cheat
if options.rom.bc_recon == 5
    basis = B;
elseif options.rom.bc_recon == 3
    phi_inhom = options.rom.phi_inhom;
    basis = [B phi_inhom];
else
    warning('not implemented')
end
V_best = basis*(basis'*(Om.*(snapshots_V_total)));
V_best2 = basis*(basis'*(Om.*(V_total)));

%% plot best approx coefficients
% if j==Nsim
%     a_best = (basis'*(Om.*(snapshots_V_total)));
%     figure(67893)
%     for k = 1:6 %size(a_best,1)
% %         plot(t_vec,a_best(10*k,:))
%         semilogy(t_vec,abs(a_best(20*k,:)))
%         hold on
%     end
% end

%%

error_V_best = V_best - snapshots_V_total;
error_V_best2 = V_best2 - snapshots_V_total;

% error_V_best_2 = weightedL2norm(error_V_best,options.grid.Om)./V_2_ref;
error_V_best_2 = weightedL2norm(error_V_best,options.grid.Om)/V_2_ref_avg;
error_V_best_22 = weightedL2norm(error_V_best,options.grid.Om)/V_2_ref_avg;

disp('velocity error')
% figure(99)
% plot(t_vec,error_V_2,color,'displayname',"ROM "+suffix);
%%
% plot(t_vec,error_V_2,'color',color,'displayname',name); %%non-botch
plot(t_vec,error_V_2,'color',color,'linestyle', linestyle, 'linewidth', linewidth,'displayname',name); %%non-botch
hold on


set(gca,'Yscale','log');
xlabel('t')
    set(gcf, 'Position', [100, 100, 400, 300])
% legend('show','NumColumns',2,'Orientation','horizontal')
legend('show','NumColumns',3,'Orientation','vertical')
    set(gcf, 'Position', [100, 100, 400, 300])
grid on



% if mod(j,2) == 0
%     plot(t_vec,error_V_2,'color',color_old,'linestyle','-' ,'displayname',name);
% else
%     plot(t_vec,error_V_2,'color',color,'linewidth',2,'linestyle',':','displayname',name);
%     color_old = color;
% end
%%
sum_error_V_2(j) = sum(error_V_2)
figure(9999) % botch
hold on
% plot(t_vec,error_V_best_2,'color',color,'linestyle','--','displayname',"best approx "+suffix);
plot(t_vec,error_V_best_2,'color',color,'linestyle','--','displayname',name);

sum_error_V_best_2(j) = sum(error_V_best_2)




% plot(t_vec,error_V_best_22,color3,'displayname',"best approx "+suffix);
% plot(t_vec,error_V_best_22,'color',color,'linestyle','--','displayname',"best approx "+suffix);
% sum_error_V_best_22(j) = sum(error_V_best_22)

%% auxiliary plots
% V_h_norm = weightedL2norm(snapshots_V_total,options.grid.Om);
% V_r_norm = weightedL2norm(V_total,options.grid.Om);
% lower_bound = abs(k - snapshots.k(snapshot_indx'))./(V_h_norm+V_r_norm);
% 
% plot(t_vec,lower_bound',color3,'displayname',"lower bound M="+M+suffix);
% plot(t_vec,abs(V_h_norm-V_r_norm),color4,'displayname',"sanity check M="+M+suffix);

%%

set(gca,'Yscale','log');
xlabel('t')
    set(gcf, 'Position', [100, 100, 400, 300])
% legend('show','NumColumns',2,'Orientation','horizontal')
legend('show','NumColumns',3,'Orientation','vertical')
    set(gcf, 'Position', [100, 100, 400, 300])
    ylim([1e-9 1])
grid on