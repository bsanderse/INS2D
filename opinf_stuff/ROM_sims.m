function [A_ROMs,A_dot_ROMs,A_raw_ROMs] = ROM_sims(D,C,a0s,dt,Nt,skip)

[M,n_trajes] = size(a0s);

A_raw_ROMs = zeros(M,Nt,n_trajes);

n_snapshots = length(1:skip:Nt-1);
A_ROMs = zeros(M,n_snapshots,n_trajes);
A_dot_ROMs = zeros(M,n_snapshots,n_trajes);
for i =1:n_trajes
    a0 = a0s(:,i);
    A_raw_ROM = ROM_sim(D,C,a0,dt,Nt);
    [A_dot_ROM,A_ROM] = time_difference_quotient(A_raw_ROM, "forward euler",dt,skip);
    A_ROMs(:,:,i) = A_ROM;
    A_dot_ROMs(:,:,i) = A_dot_ROM;
    A_raw_ROMs(:,:,i) = A_raw_ROM;
end

A_ROMs = A_ROMs(:,:);
A_dot_ROMs = A_dot_ROMs(:,:);


