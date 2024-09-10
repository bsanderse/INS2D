function error = relative_state_error_best_approx(FOM_snapshots,ROM_basis,options)
% relative state error as defined in 
% https://arxiv.org/pdf/2401.02889.pdf#cite.kaptanoglu2021promoting
% equation (18)

error = norm(FOM_snapshots - ROM_basis*(ROM_basis'*(options.grid.Om.*FOM_snapshots)),"fro")/norm(FOM_snapshots,"fro");