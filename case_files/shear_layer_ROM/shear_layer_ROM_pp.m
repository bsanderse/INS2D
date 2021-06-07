%% post-processing shear layer
Npx = options.grid.Npx;
Npy = options.grid.Npy;

% %% kinetic energy
% figure
% plot(time,(k-k(1))/k(1),'s-')
% grid
% xlabel('t');
% ylabel('(k(t)-k(0))/k(0)');
% 
% %% vorticity
% % compare e.g. with PhD thesis figure 3.4
% 
% figure
% omega = get_vorticity(V,t,options);
% omega = reshape(omega,Npx+1,Npy+1);
% % for Re=1000: labels = -4:0.5:4;
% labels= 20;
% contour(x,y,omega',labels);
% axis square
% colorbar
% grid

colors = ['b' 'r' 'g' 'k' 'b' 'r' 'g' 'k'];
% lines = ['-' '-' '-' '-' '--' '--' '--' '--'];
% 
% lines = ['-';'--'];
lines = '--';

%% Fig. 3. Singular values for inviscid shear-layer roll-up
% if j==1
%     singulars = figure;
% end
% if j==1 || j==2
%     if (size(S,2)>1)
%         Sigma = diag(S);
%     else
%         Sigma = S;
%     end
%     if j==1
%         singulars = figure;
%         Sigma0 = Sigma;
%         dispname = "standard SVD";
%     else
%         dispname = "momentum-conserving SVD";
%         hold on
%     end
%     Sigma = Sigma(1:80);
%     figure(singulars);
%     semilogy(Sigma/Sigma(1),'s','displayname', dispname);
%     hold on
%     title('singular values')
%     legend('show')
% end


% [~,S0] = getBasis(V_svd,options,M);
% if (size(S0,2)>1)
%     Sigma0 = diag(S0);
% else
%     Sigma0 = S0;
% end
% semilogy(Sigma0/Sigma0(1),'s');

%% Fig. 4 
k_r = k;
umom_r = umom;
% load("../../"+snapshot_data,'k','umom')
if (~exist('test_data'))
    test_data = snapshot_data;
end
load("../../"+test_data,'k','umom')


if j==1
    energy1 = figure;
    % title('energy 1')
    % hold on
    energy2 = figure;
    % title('energy 2')
    % hold on
    momentum = figure;
    % title('momentum')
    % hold on
end
% if j>4
% %     suffix = " explicit";
%     suffix = " mc";
%     lines = "--";
% else
% %     suffix = " implicit";
%     suffix = " ";
%     lines = "-";
% end

figure(energy1)
semilogy(abs(k_r-k(1))/k(1),'color', colors(j), 'LineStyle', lines, ...
    'displayname',"M="+M+suffix)
hold on
title('energy error')
xlabel('t')
ylabel("(K_r^n-K_h(0))/K_h(0)")
legend('show')

figure(momentum)
semilogy(abs(umom_r-umom(1))/umom(1),'color', colors(j), ...
    'LineStyle', lines,'displayname',"M="+M+suffix)
hold on
title('momentum error')
xlabel('t')
ylabel('error in u component of global momentum')
legend('show')

figure(energy2)
semilogy(abs(k_r-k_r(1))/k_r(1),'color', colors(j), ...
    'LineStyle', lines, 'displayname',"M="+M+suffix)
hold on
title('energy error')
xlabel('t')
ylabel("(K_r^n-K_r^0)/K_r^0")
legend('show')
