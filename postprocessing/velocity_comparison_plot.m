V_total = [uh_total vh_total]'; %[uh_total; vh_total];
snapshots_V_total = [test_snapshots.uh_total(snapshot_indx,:) test_snapshots.vh_total(snapshot_indx,:)]';

%% wrong
% best possible approximation given the projection:
% note that the norm should be consistent with the
% optimization problem used in the SVD

% B = options.rom.B;
% Om = options.grid.Om;
% if options.verbosity.equivalence_cheat
%     V_best = B*(B'*(Om.*(snapshots_V_total)));
% else
%     phi_inhom = options.rom.phi_inhom;
% %     V_best = B*(B'*(Om.*(snapshots_V_total)))+snapshots.Vbc; %actually I probably should also take the best approximation of Vbc onto phi_inhom here
%     V_best = B*(B'*(Om.*(snapshots_V_total))) ...
%             +phi_inhom*phi_inhom'*(Om.*snapshots.Vbc); 
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

%%

if mod(j,2) == 1
    V_total_old = V_total;
    V_best_old = V_best;
    suffix_old = suffix;
    M_old = M;
    name_old = name;
    color_old = color;
%     color2_old = color2;
else
    error_V = V_total - V_total_old;
    % inf-norm
    error_V_inf = max(abs(error_V),[],1);
    
    V_ref    = ones(size(snapshots_V_total));
    
    V_2_ref  = weightedL2norm(V_ref,options.grid.Om);
    % error_V_2 = weightedL2norm(error_V,options.grid.Om)./V_2_ref;
    error_V_2 = weightedL2norm(error_V,options.grid.Om);

    error_V_best = V_best - V_best_old;
%     errors_V_best(j) = norm(error_V_best)

    % error_V_best_2 = weightedL2norm(error_V_best,options.grid.Om)./V_2_ref;
    error_V_best_2 = weightedL2norm(error_V_best,options.grid.Om);

    disp('rom equivalence')
    figure(999)
%     plot(t_vec,error_V_2,'color',color_old,'displayname',"R = "+num2str(M_old));
%     plot(t_vec,error_V_2,'color',color_old,'displayname',name);
    plot(t_vec,error_V_2,'displayname',name + "/" + name_old);
%     plot(t_vec,error_V_2,color_old,'displayname',"velo diff: "+ name +" / "+ name_old);
    sum2_error_V_2(j) = sum(error_V_2)
    hold on
%     plot(t_vec,error_V_best_2,color2_old,'displayname',"best approx diff: "+ name +" / "+ name_old);
    sum2_error_V_best_2(j) = sum(error_V_best_2)
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
ylabel("||\cdot||_{\Omega_h}")
% legend('show','NumColumns',2,'Orientation','horizontal')
legend('show','NumColumns',1,'Orientation','horizontal')

% legend('velocity difference mthesis ROMs Mhom=2 Minhom=2 Mbc=2 Mp=2', ...
% 'velocity best approx difference mthesis ROMs Mhom=2 Minhom=2 Mbc=2 Mp=2', ...
% 'velocity difference optimal ROMs Mhom=2 Minhom=2 Mbc=2 Mp=2', ...
% 'velocity best approx difference optimal ROMs Mhom=2 Minhom=2 Mbc=2 Mp=2', ...
% 'velocity difference qr ROMs Mhom=2 Minhom=2 Mbc=2 Mp=2', ...
% 'velocity best approx difference qr ROMs Mhom=2 Minhom=2 Mbc=2 Mp=2', ...
% 'NumColumns',1,'Orientation','horizontal')