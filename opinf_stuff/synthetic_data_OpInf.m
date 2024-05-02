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
K = 2;
Uhats = zeros(r+r^2,r,r,K);
RHSs = zeros(r,r,r,K);
for k=1:K
    for i=1:r
        U = zeros(r,1);
        U(i) = k;
        for j =1:r
            V = zeros(r,1);
            V(j) = k;
            VU = kron(V,U);

            Uhats(:,i,j,k) = [U; VU];
            RHSs(:,i,j,k) = nu*D*U - C*VU;
        end
    end
end
Uhats = Uhats(:,:);
RHSs = RHSs(:,:);

%%

[nuD,C_OpInf] = OpInf_core(Uhats,RHSs);


D_OpInf = nuD/nu;