clear all;
close all;
% toy problems for op inf - 1D Burgers periodic BC
% function D_error = opinf_test_burgers(nt)

dx = .1;
N = 3;

D = -2*diag(ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
D(N,1) = 1;
D(1,N) = 1;
D = D/dx^2;

D = 0*D;

index = @(x) mod(x-1,N)+1;
kron_ind = @(i,j) i+(j-1)*N;

C = zeros(N,N^2);
for i = 1:N
    C(i,kron_ind(index(i-1),index(i-1))) = 1;
    C(i,kron_ind(index(i+1),index(i+1))) = -1;
    C(i,kron_ind(index(i-1),i)) = 1;
    C(i,kron_ind(index(i+1),i)) = -1;
end

C = C/(3*dx);

% C = 0*C;
 
nu = .01;
dt = .0001;

%% start simulation
n_offsets = 1;
offsets = 1:n_offsets;
n_x0 = N;
[V,E] = eig(D);
X0 = V + permute(repmat(offsets,N,1,n_x0),[3 1 2]);
X0 = X0(:,:);
n_x0 = n_x0*n_offsets;

% n_x0 = 1;
% % x0 = (sin(2*pi*(1:N)/N)').^2;
% X0 = (sin(2*pi*(1:N)/N)') + 2;

% v_max = max(X0,[],'all');
% 
% dt = .1*dx/v_max;
% dt = dx/v_max;


nts = 10*(1:15);

for tt = 1:numel(nts)
nt = nts(tt);
% nt = 40; % number of time steps to be performed per initial condition minus 1 (for initial condition)
           % = number of snapshots per trajectory

K = n_x0 * nt;

X = zeros(N,nt,n_x0);
X_dot = zeros(N,nt-1,n_x0);

energies = zeros(nt,n_x0);
% X(:,1) = x0;

for ii = 1:n_x0
    x = X0(:,ii);
    X(:,1,ii) = x;
    energies(1,ii) = norm(x)^2/2;
    for i = 1:nt-1
        %% forward Euler
        x = x + dt*(nu*D*x + C*kron(x,x));  
        %% Gaus method (energy-conserving!)
            % still waiting for license
        %% linear implicit method
%             T = (nu*D + C*kron(x,eye(N)));
%             x12 = (eye(N) - .5*dt*T)\x;
%             x = x + dt*T*x12;
% %             X_dot(:,i,ii);
        %%
        X(:,1+i,ii) = x;
        energies(1+i,ii) = norm(x)^2/2;
    end
end

%% end simulation

[Diff,Conv] = standardOpInf(X,dt,nu);
Conv_r = reduced_convection_operator(Conv);
C_r = reduced_convection_operator(C);
[Diff2,Conv2] = standardOpInf2(X,dt,nu);

D_error(tt) = norm(D-Diff)
D_error2(tt) =norm(D-Diff2)

C_error(tt) = norm(reduced_convection_operator(C) - reduced_convection_operator(Conv))
C_error2(tt) =norm(reduced_convection_operator(C) - reduced_convection_operator(Conv2))

end

figure
heatmap(X(:,:))
grid off
% 
% figur ("show")

figure
semilogy(nts,D_error,'gx-',displayname = "||D\_intrusive - D\_OpInf||")
hold on
semilogy(nts,D_error2,'gx:',displayname = "||D\_intrusive - D\_OpInf2||")

% title("||D\_intrusive - D\_OpInf||")
% xlabel("number of snapshots per initial condition trajectory")
% ylabel("error magnitude")

% figure
semilogy(nts,C_error,'bx-',displayname = "||C\_intrusive - C\_OpInf||")
semilogy(nts,C_error2,'bx:',displayname = "||C\_intrusive - C\_OpInf2||")
% title("||C\_intrusive - C\_OpInf||")
xlabel("number of snapshots per initial condition trajectory")
ylabel("error magnitude")
legend("show")

