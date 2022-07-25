V_total = [uh_total vh_total]'; %[uh_total; vh_total];
snapshots_V_total = [test_snapshots.uh_total(snapshot_indx,:) test_snapshots.vh_total(snapshot_indx,:)]';
error_V = V_total - snapshots_V_total;
% inf-norm
error_V_inf = max(abs(error_V),[],1);

V_ref    = ones(size(snapshots_V_total));

V_2_ref  = weightedL2norm(V_ref,options.grid.Om);
% error_V_2 = weightedL2norm(error_V,options.grid.Om)./V_2_ref;
error_V_2 = weightedL2norm(error_V,options.grid.Om);

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
if options.verbosity.equivalence_cheat
    basis = B;
else
    phi_inhom = options.rom.phi_inhom;
    basis = [B phi_inhom];
end
V_best = basis*(basis'*(Om.*(snapshots_V_total)));

%%

error_V_best = V_best - snapshots_V_total;

% error_V_best_2 = weightedL2norm(error_V_best,options.grid.Om)./V_2_ref;
error_V_best_2 = weightedL2norm(error_V_best,options.grid.Om);

L2_V_best_errors(j) = sum(error_V_best_2)

figure(99)
plot(t_vec,error_V_2,color,'displayname',"ROM M="+M+suffix); 
hold on
plot(t_vec,error_V_best_2,color2,'displayname',"best approx M="+M+suffix);

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
ylabel('velocity error')
legend('show','NumColumns',2,'Orientation','horizontal')