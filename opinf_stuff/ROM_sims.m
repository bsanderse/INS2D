function [A_ROMs,A_dot_ROMs,A_raw_ROMs] = ROM_sims(D,C,a0s,dt,Nt)

[M,n_trajes] = size(a0s);

A_ROMs = zeros(M,Nt-1,n_trajes);
A_dot_ROMs = zeros(M,Nt-1,n_trajes);
A_raw_ROMs = zeros(M,Nt,n_trajes);
for i =1:n_trajes
    a0 = a0s(:,i);
    A_raw_ROM = ROM_sim(D,C,a0,dt,Nt);
    [A_dot_ROM,A_ROM] = time_difference_quotient(A_raw_ROM, "forward euler",dt);
    A_ROMs(:,:,i) = A_ROM;
    A_dot_ROMs(:,:,i) = A_dot_ROM;
    A_raw_ROMs(:,:,i) = A_raw_ROM;
end

A_ROMs = A_ROMs(:,:);
A_dot_ROMs = A_dot_ROMs(:,:);


