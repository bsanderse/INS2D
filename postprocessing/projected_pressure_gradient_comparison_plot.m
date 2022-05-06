G_h = -options.discretization.M';

snapshots_p_total = test_snapshots.p_total(snapshot_indx,:)';
G_p = G_h*p_total';

    if options.rom.bc_recon == 3
        phi = [options.rom.B options.rom.phi_inhom];
        phi_hom = options.rom.B ;
        phi_inhom = options.rom.phi_inhom;
    elseif options.rom.bc_recon == 5
        phi = options.rom.B;
    end

if mod(j,2)==1
    G_p_old = G_p;
    name_old = name;
else
    error_G_p = G_p-G_p_old;

    if options.rom.bc_recon == 3
        phi = [options.rom.B options.rom.phi_inhom];
    elseif options.rom.bc_recon == 5
        phi = options.rom.B;
    end

    error_phi_G_p = phi'*error_G_p;

    M = options.rom.M;
    one_vec = ones(M,1);

    error_phi_G_p_2 = weightedL2norm(error_phi_G_p, one_vec);

    %%
        error_phi_hom_G_p = phi_hom'*error_G_p;
        error_phi_inhom_G_p = phi_inhom'*error_G_p;

    one_vec1 = ones(size(phi_hom,2),1);
    one_vec2 = ones(size(phi_inhom,2),1);

    error_phi_hom_G_p_2 = weightedL2norm(error_phi_hom_G_p, one_vec1);
    error_phi_inhom_G_p_2 = weightedL2norm(error_phi_inhom_G_p, one_vec2);
    

    %%

    figure(9812)
    plot(t_vec,error_phi_G_p_2,color,'displayname',"pres diff: "+ name +" / "+ name_old);

    hold on

end

set(gca,'Yscale','log');
xlabel('t')
ylabel('projected pressure gradient difference')
% legend('show','NumColumns',2,'Orientation','horizontal')
% legend('show')