snapshots_p_total = test_snapshots.p_total(snapshot_indx,:)';
error_p = p_total' - snapshots_p_total;

Np = options.grid.Np;
one_vec = ones(Np,1);

error_p_2 = weightedL2norm(error_p, one_vec);

figure(98)
plot(t_vec,error_p_2,color,'displayname',"ROM M="+M+suffix);
hold on

if options.rom.bc_recon == 5
% best possible approximation given the projection:
% note that the norm should be consistent with the
% optimization problem used in the SVD
Bp = options.rom.Bp;
Omp         = options.grid.Omp;
p_best = Bp*(Bp'*(Omp.*snapshots_p_total));

error_p_best = p_best - snapshots_p_total;

error_p_best_2 = weightedL2norm(error_p_best, one_vec);
 
plot(t_vec,error_p_best_2,color2,'displayname',"best approx M="+M+suffix);
end

set(gca,'Yscale','log');
xlabel('t')
ylabel('pressure error')
legend('show','NumColumns',2,'Orientation','horizontal')
% legend('show')