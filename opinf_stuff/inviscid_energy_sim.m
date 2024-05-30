function [energies,X] = inviscid_energy_sim(C,dt,x0,nt,method)

x = x0;
N = numel(x0);

X = zeros(N,nt);
X(:,1) = x;
energies = zeros(nt,1);
energies(1) = norm(x)^2/2;

for i = 1:nt-1
    switch method
        case "linear implicit"
            T = C*kron(x,eye(N));
            x12 = (eye(N) - .5*dt*T)\x;
            x = x + dt*T*x12;
        case "Gauss method"
            options = optimset('Display','off');
            x = fsolve(@(x1) x + dt*C*kron(x1,x1) - x1,x,options);
    end    

    % %% forward Euler
    % x = x + dt*(nu*D*x + C*kron(x,x));
    % %% Gaus method (energy-conserving!)
    % % options = optimset('Display','off');
    % % x = fsolve(@(x1) x + dt*(nu*D*x1 + C*kron(x1,x1)) - x1,x,options);
    % %% linear implicit method
    % %             T = (nu*D + C*kron(x,eye(N)));
    % %             x12 = (eye(N) - .5*dt*T)\x;
    % %             x = x + dt*T*x12;
    % % %             X_dot(:,i,ii);
    %%
    X(:,1+i) = x;
    energies(1+i) = norm(x)^2/2;
end