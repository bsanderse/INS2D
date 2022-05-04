V_total = [uh_total vh_total]'; %[uh_total; vh_total];
snapshots_V_total = [test_snapshots.uh_total(snapshot_indx,:) test_snapshots.vh_total(snapshot_indx,:)]';


% best possible approximation given the projection:
% note that the norm should be consistent with the
% optimization problem used in the SVD
B = options.rom.B;
Om = options.grid.Om;
if options.verbosity.equivalence_cheat
    V_best = B*(B'*(Om.*(snapshots_V_total)));
else
    V_best = B*(B'*(Om.*(snapshots_V_total)))+snapshots.Vbc;
end

if mod(j,2) == 1
    V_total_old = V_total;
    V_best_old = V_best;
    suffix_old = suffix;
    M_old = M;
else
    error_V = V_total - V_total_old;
    % inf-norm
    error_V_inf = max(abs(error_V),[],1);
    
    V_ref    = ones(size(snapshots_V_total));
    
    V_2_ref  = weightedL2norm(V_ref,options.grid.Om);
    % error_V_2 = weightedL2norm(error_V,options.grid.Om)./V_2_ref;
    error_V_2 = weightedL2norm(error_V,options.grid.Om);

    error_V_best = V_best - V_best_old;

    % error_V_best_2 = weightedL2norm(error_V_best,options.grid.Om)./V_2_ref;
    error_V_best_2 = weightedL2norm(error_V_best,options.grid.Om);

    figure(999)
    plot(t_vec,error_V_2,color,'displayname',"deviation ROM M="+M+suffix +"/ ROM M=" +M_old+suffix_old);
    hold on
    plot(t_vec,error_V_best_2,color2,'displayname',"deviation best approx M="+M+suffix +"/ best approx M=" +M_old+suffix_old);

end

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