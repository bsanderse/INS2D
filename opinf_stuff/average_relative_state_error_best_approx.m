function error = average_relative_state_error_best_approx(Nt,FOM_snapshots_s,ROM_basis,options)
% average relative state error as defined in 
% https://arxiv.org/pdf/2401.02889.pdf#cite.kaptanoglu2021promoting
% equation (18)
% assuming Forward Euler time discretization!!!


[NV,Nt_n_trajes] = size(FOM_snapshots_s);
n_trajes = Nt_n_trajes / Nt;
FOM_snapshots_tensor = reshape(FOM_snapshots_s,NV,Nt,n_trajes);

error = 0;
for i=1:n_trajes
    FOM_snapshots = FOM_snapshots_tensor(:,:,i);

    error = error + relative_state_error_best_approx(FOM_snapshots,ROM_basis,options);
end

error = error/n_trajes;
