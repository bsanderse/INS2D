%% ugliest botch ever
figure(3768)
a = plot([2 3]);
hold on
b = plot(2*[2 3]);
c = plot(3*[2 3]);
d = plot(4*[2 3]);
e = plot(5*[2 3]);
f = plot(6*[2 3]);
g = plot(7*[2 3]);

colors = [a.Color; b.Color; c.Color; d.Color; e.Color; f.Color; g.Color];
% colors = kron(colors,[1; 1]);

% color = colors(j,:);

%%

ts = t_start + [0:nt]*dt;

err1 = zeros(nt+1,1);
err2 = zeros(nt+1,1);
ybc_norms = zeros(nt+1,1);

% Rs = [2 5 10 20 40 79 80];
% Rs = [2 5 10 20 40]; %79 80];
% Rs = 20 + [1:6];

% Rs = [2 5 10 20 23 24 40];
% Rs = 80;

Rs = [2 5 10 20 22 23 24];
% Rs = [10 20 40 80 81 85];
% Rs = [81 82 83 84 85];
% Rs = [2 5 20 79 80 84];


phi_bc0 = options.rom.phi_bc;

for l = 1:numel(Rs)
    R = Rs(l);

    for k = 1:nt+1
        t = ts(k);
        y_bc0 = get_bc_vector_yBC(t,options);
        phi_bc = phi_bc0(:,1:R);
        BC_DEIM = options.rom.BC_DEIM;
        options.rom.BC_DEIM = 0;
        a_bc1 = get_a_bc(t,options);
        y_bc1 = phi_bc*a_bc1(1:R);
        options.rom.BC_DEIM = 1;
        options.rom.BC_DEIMdim = R;
        [options.rom.bc_deim_inds, options.rom.bc_deim_PTUinv] = getDEIMinds(phi_bc,options.rom.BC_DEIMdim);
        options.rom.ybc_coords = ybc_coords(options.rom.bc_deim_inds,t,options);
        a_bc2 = get_a_bc(t,options);
        y_bc2 = phi_bc*a_bc2;
        options.rom.BC_DEIM = BC_DEIM;

        err1(k) = norm(y_bc0-y_bc1);
        err2(k) = norm(y_bc0-y_bc2);
        ybc_norms(k) = norm(y_bc0);
%         figure(467)
%         plot(y_bc2-y_bc0)
%         hold on
    end

    ybc_norm_avg = sum(ybc_norms)/(nt+1);

    linewidth = 1;
    figure(377111)
    color = colors(l,:);
    semilogy(ts,err1/ybc_norm_avg,'color',color,'linewidth',linewidth,'linestyle','--', 'displayname', "R="+num2str(R))
    hold on
    semilogy(ts,err2/ybc_norm_avg,'color',color,'linewidth',linewidth,'linestyle','-', 'displayname', "R="+num2str(R))


end

% legend('show')
    set(gcf, 'Position', [100, 100, 400, 600])
    grid on
legend('show','NumColumns',3,'Orientation','vertical')