%% post-processing shear layer
Npx = options.grid.Npx;
Npy = options.grid.Npy;

%% kinetic energy
figure(13)
% plot(time,(k-k(1))/k(1),'s-','displayname',"dt = "+num2str(dt))
plot(time,(k-k(1))/k(1),'s-')
% semilogy(time,abs((k-k(1))/k(1)),'s-','displayname',"dt = "+num2str(dt))

hold on
grid
xlabel('t');
ylabel('(k(t)-k(0))/k(0)');
legend('show')

%% kinetic energy error
figure(14)
% plot(time,(k-k(1))/k(1),'s-','displayname',"dt = "+num2str(dt))
plot(time,abs(k-k(1))/k(1),'s-')
% semilogy(time,abs((k-k(1))/k(1)),'s-','displayname',"dt = "+num2str(dt))

hold on
grid
xlabel('t');
ylabel('abs(k(t)-k(0))/k(0)');
legend('show')

%% vorticity
% compare e.g. with PhD thesis figure 3.4

% figure
% omega = get_vorticity(V,t,options);
% omega = reshape(omega,Npx+1,Npy+1);
% % for Re=1000: labels = -4:0.5:4;
% labels= 20;
% contour(x,y,omega',labels);
% axis square
% colorbar
% grid