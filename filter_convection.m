% construct filter for convective terms
function u_filtered = filter_convection(u,diff_matrix,bc,alfa)

    u_filtered = u + alfa*(diff_matrix*u + bc);


end