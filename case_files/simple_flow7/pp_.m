% disp('pp.m is empty')

%% kinetic energy
figure
plot(t_start:dt:t_end,(k-k(1))/k(1),'s-');
title('normalized kinetic energy evolution')
if save_file == 1
    savefig(strcat('../../',path_results,'/ kinetic energy.fig'))
end


%% divergence of velocity field
figure
semilogy(t_start:dt:t_end,maxdiv);
title('divergence of velocity field')
if save_file == 1
    savefig(strcat('../../',path_results,'/ divergence.fig'))
end

