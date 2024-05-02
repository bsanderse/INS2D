% try operator inference in INS2D framework but with synthetic data
function [D_OpInf,C_OpInf] = synthetic_FOM_data_OpInf(nu,D,C,phi)

[N,r] = size(phi);

%% 
% Us = zeros(r);
% RHSs = zeros(r);
% for i=1:r
%     U = zeros(r,1);
%     U(i) = 1;
%     Us(:,i) = U;
%     RHSs(:,i) = nu*D*U - C*kron(U,U);
% end

%%

sparse_ = true;
K = 2;
% if sparse_
%     Uhats = sparse(r+r^2,r*r*K);
%     RHSs = sparse(r,r*r*K);
% else
    Uhats = zeros(r+r^2,r,r,K);
    RHSs = zeros(r,r,r,K);
% end

for k=1:K
    for i=1:r
        U = sparse(i,1,k,N,1);
        for j =1:r
            V = sparse(j,1,k,N,1);
            %% heterogenous data
            % VU = kron(V,U); 
            % Uhats(:,i,j,k) = [U; VU];
            % RHSs(:,i,j,k) = nu*D*U - C*VU;
            %% homogeneous data
            U = U + V;
            UU = kron(U,U);
            a = phi'*U;
            Uhats(:,i,j,k) = [a;kron(a,a)];
            RHSs(:,i,j,k) = phi'*(nu*D*U - C*UU);
            %%
        end
    end
end

    % if ~sparse_
        Uhats = Uhats(:,:);
        RHSs = RHSs(:,:);
    % end

%%

[nuD,C_OpInf] = OpInf_core(Uhats,RHSs);
D_OpInf = nuD/nu;
