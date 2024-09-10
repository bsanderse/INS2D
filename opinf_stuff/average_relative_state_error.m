function error = average_relative_state_error(D,C,a0s,dt,Nt,FOM_snapshots_s,ROM_basis,skip)
% average relative state error as defined in 
% https://arxiv.org/pdf/2401.02889.pdf#cite.kaptanoglu2021promoting
% equation (18)
% assuming Forward Euler time discretization!!!

[~,~,A_ROMs_s] = ROM_sims(D,C,a0s,dt,Nt,skip);

[r,n_trajes] = size(a0s);

A_ROMs_tensor = reshape(A_ROMs_s,r,Nt,n_trajes);

NV = size(FOM_snapshots_s,1);
FOM_snapshots_tensor = reshape(FOM_snapshots_s,NV,Nt,n_trajes);

error = 0;
for i=1:n_trajes
    ROM_snapshots = A_ROMs_tensor(:,:,i);
    FOM_snapshots = FOM_snapshots_tensor(:,:,i);

    error = error + relative_state_error(FOM_snapshots,ROM_snapshots,ROM_basis);
end

error = error/n_trajes;
