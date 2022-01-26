% Jacobian check
% [maxres,Fres,dF] = F(V,C,p,t,options,getJacobian,nopressure)

%%
% func = @(V) F(V,V,p,t,options,1,0);

% N = 3;
% Conv_diag = magic(3).^2;
% V = [3 7 -2]';
% V0 = [3 7 -2]';


% Conv_diag = options.grid.C;
% func = @(V) (Conv_diag*V).*V; % correct one
% func = @(V) (Conv_diag*V0).*V;
% func = @(V) (Conv_diag*V).*V0;
% func = @(V) (Conv_diag*V).*V0 + (Conv_diag*V0).*V;

% Jac = diag(Conv_diag*V) + diag(V)*Conv_diag; % correct one
% Jac = diag(Conv_diag*V0) + diag(V0)*Conv_diag;
% Jac = Conv_diag.*V0;
% Jac = diag(V0)*Conv_diag;
% Jac = diag(Conv_diag*V0);
%%
% id_normal = options.grid.id_normal;
% id_tangential = options.grid.id_tangential;
% id_n_t = id_normal+id_tangential;
% % V_n_t = id_n_t.*V;
% gO = options.BC.gO;
% dgO = options.BC.dgO;
% 
% % func = @(V) - (id_n_t.*V).*gO(V);
% % Jac =  - diag(gO(V).*id_n_t + (id_n_t.*V).*dgO(V));
% 
% % func = @(V) - diag(id_n_t)*diag(gO(V))*V;
% % Jac = - diag(id_n_t)*diag(gO(V)) - diag(id_n_t)*diag(dgO(V))*diag(V);
% 
% V0 = V;
% % func = @(V) - diag(id_n_t)*diag(gO(V0))*V;
% % Jac = - diag(id_n_t)*diag(gO(V0));
% 
% func = @(V) - diag(id_n_t)*diag(gO(V))*V0;
% Jac = - diag(id_n_t)*diag(dgO(V))*diag(V0);


%%

[~,F0,Jac] = F(V,V,p,t,options,1,0);
% F0 = func(V);



epsilons = 10.^(-1*(-4:8));
% epsilons = 10.^(-3);

errors = 0*epsilons;
for jj = 1:numel(epsilons)
epsilon = epsilons(jj);
Jac2 = zeros(numel(V));

for i = 1:numel(V)
    Vi = V;
    Vi(i) = Vi(i)+epsilon;
    
    [~,Fi,~] = F(Vi,Vi,p,t,options,1,0);
%     Fi = func(Vi);
    
    dFi = (Fi-F0)/epsilon;
    Jac2(i,:) = dFi;
end

errorj = norm(Jac-Jac2');
errors(jj) = errorj;

end
figure
loglog(epsilons,errors,'x')
1



