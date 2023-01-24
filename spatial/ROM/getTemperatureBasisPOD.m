function [BT,MT] = getTemperatureBasisPOD(snapshots,sample_index,options)
% build basis for velocity field based on snapshot data

disp('computing SVD of pressure snapshots...');
svd_start2 = toc;
% note p_total is stored as a Nt*Np matrix instead of Np*Nt which we use for
% the SVD
% use same snapshot_indx that was determined for velocity

% select snapshots
T_total_snapshots = snapshots.T_total';

T_svd  = T_total_snapshots(:,sample_index);

% take first Mp columns of Wp as a reduced basis
% (better is to look at the decay of the singular values in Sp to determine M)
if (isfield(options.rom,'MT'))
    MT = options.rom.MT;
else
    % if not defined, use same number of modes as for velocity
    warning('number of temperature modes not defined, defaulting to number of velocity modes');
    MT = options.rom.MT;
end

if (options.rom.weighted_norm_T == 0)

    [WT,ST] = getBasis(T_svd,options,options.rom.MT);

    % perform SVD
    %     [Wp,Sp,Zp] = svd(p_svd,'econ');

elseif (options.rom.weighted_norm_T == 1)

    Np          = options.grid.Np;
    Omp         = options.grid.Omp;
    Omp_sqrt    = spdiags(sqrt(Omp),0,Np,Np);
    Omp_invsqrt = spdiags(1./sqrt(Omp),0,Np,Np);

    % make weighted snapshot matrix
    Tmod = Omp_sqrt*T_svd;

    % getBasis can use different methods to get basis: SVD/direct/snapshot
    % method
    [WT,ST] = getBasis(Tmod,options);

    % transform back
    WT = Omp_invsqrt*WT;
end

BT = WT(:,1:MT);
options.rom.BT = BT;

svd_end = toc - svd_start2

figure(21)
hold on
if (size(ST,2)>1)
    SigmaT = diag(ST);
else
    SigmaT = ST;
end
semilogy(SigmaT/SigmaT(1),'o');

