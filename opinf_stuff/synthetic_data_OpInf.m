% try operator inference in INS2D framework but with synthetic data
function [D_OpInf,C_OpInf] = synthetic_data_OpInf(nu,D,C)

r = size(D,1);

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

sparse_ = false;
K = 2;
if sparse_
    Uhats = sparse(r+r^2,r*r*K);
    RHSs = sparse(r,r*r*K);
else
    Uhats = zeros(r+r^2,r,r,K);
    RHSs = zeros(r,r,r,K);
end

for k=1:K
    for i=1:r
        U = zeros(r,1);
        U(i) = k;
        for j =1:r
            V = zeros(r,1);
            V(j) = k;
            %% heterogenous data
            % VU = kron(V,U); 
            % Uhats(:,i,j,k) = [U; VU];
            % RHSs(:,i,j,k) = nu*D*U - C*VU;
            %% homogeneous data
            U = U + V;
            UU = kron(U,U);
            if sparse_ 
                Uhats(:,i+(j-1)*r+(k-1)*K) = [U; UU];
                RHSs(:,i+(j-1)*r+(k-1)*K) = nu*D*U - C*UU;
            else
                Uhats(:,i,j,k) = [U; UU];
                RHSs(:,i,j,k) = nu*D*U - C*UU;
            end
            %%
        end
    end
    end
    if ~sparse_
        Uhats = Uhats(:,:);
        RHSs = RHSs(:,:);
    end

%%

[nuD,C_OpInf] = OpInf_core(Uhats,RHSs);
D_OpInf = nuD/nu;