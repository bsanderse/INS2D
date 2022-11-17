function [Bp,Mp] = getPressureBasisPOD(snapshots,sample_index,options)
% build basis for velocity field based on snapshot data

disp('computing SVD of pressure snapshots...');
svd_start2 = toc;
% note p_total is stored as a Nt*Np matrix instead of Np*Nt which we use for
% the SVD
% use same snapshot_indx that was determined for velocity

% select snapshots
p_total_snapshots = snapshots.p_total';
if (options.rom.pressure_mean == 1)
    % subtract temporal mean
    options.rom.p_mean = mean(p_total_snapshots,2);
    p_total_snapshots = p_total_snapshots - options.rom.p_mean;
end
p_svd  = p_total_snapshots(:,sample_index);

% take first Mp columns of Wp as a reduced basis
% (better is to look at the decay of the singular values in Sp to determine M)
if (isfield(options.rom,'Mp'))
    Mp = options.rom.Mp;
else
    % if not defined, use same number of modes as for velocity
    warning('number of pressure modes not defined, defaulting to number of velocity modes');
    Mp = options.rom.M;
end

if (options.rom.weighted_norm == 0)

    [Wp,Sp] = getBasis(p_svd,options,options.rom.Mp);

    % perform SVD
    %     [Wp,Sp,Zp] = svd(p_svd,'econ');

elseif (options.rom.weighted_norm == 1)

    Np          = options.grid.Np;
    Omp         = options.grid.Omp;
    Omp_sqrt    = spdiags(sqrt(Omp),0,Np,Np);
    Omp_invsqrt = spdiags(1./sqrt(Omp),0,Np,Np);

    % make weighted snapshot matrix
    pmod = Omp_sqrt*p_svd;

    % getBasis can use different methods to get basis: SVD/direct/snapshot
    % method
    [Wp,Sp] = getBasis(pmod,options);

    % transform back
    Wp = Omp_invsqrt*Wp;
end

Bp = Wp(:,1:Mp);
options.rom.Bp = Bp;

svd_end = toc - svd_start2

figure(21)
hold on
if (size(Sp,2)>1)
    SigmaP = diag(Sp);
else
    SigmaP = Sp;
end
semilogy(SigmaP/SigmaP(1),'o');

