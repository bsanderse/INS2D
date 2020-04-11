function options = operator_rom(Vbc,options)

B  = options.rom.B;
NV = options.grid.Nu+options.grid.Nv;
Z  = zeros(NV,1);

if (options.rom.weighted_norm == 0)
    Diag = options.grid.Om_inv;
elseif (options.rom.weighted_norm == 1)
    Diag = ones(NV,1);
end

% we construct the ROM operators in the following 'non-intrusive' way, i.e.
% we hardly require any knowledge of the exact FOM discretization. we
% simply make calls to diffusion and convection with the modes phi as input
% arguments

% for diffusion:
% the ROM discretization will read Diff*R + yDiff, where Diff is the reduced
% diffusion operator
% we get this by multiplying B'*Diffusion(B*R), where a call to Diffusion(V) returns
% D*V + yD
% this also works for the case that we have
% B'*Diffusion(B*R + Vbc) = 
% B'*(D*(B*R+Vbc) + yD) = (B'*D*B)*R + B'*D*Vbc + B'*yD =
% (B'*D*B)*R + B'*yDiff, where yDiff = B'*(D*Vbc + yD)
M    = options.rom.M;
Diff = zeros(M,M);

% first evaluate diffusion with zero velocity (or Vbc field) to only get boundary
% condition contribution:
% diffusion(Vbc) returns D*Vbc + yD
[d2u_bc,d2v_bc] = diffusion(Vbc,0,options,0);
yDiff           = B'*(Diag.*[d2u_bc;d2v_bc]);

% now evaluate each column of the reduced matrix as B'*D*B, where the
% boundary conditions are subtracted in order to get only the terms
% that involve the matrix
% we evaluate diffusion(B) = D*B + yD,
% and in order to get B'*D*B, we take
% diffusion(B) - diffusion(0)
[d2u_0,d2v_0] = diffusion(Z,0,options,0);
yDiff0        = B'*(Diag.*[d2u_0;d2v_0]);

for i=1:M
    [d2u,d2v] = diffusion(B(:,i),0,options,0);
    Diff(:,i) = B'*(Diag.*[d2u;d2v]) - yDiff0;
end

options.rom.Diff  = Diff;
options.rom.yDiff = yDiff;

% the convective terms are more involved
% the convective terms are of the following form:
% Cx * ( (I*V + yI).*(A*V + yA) ), where Cx is a differencing
% operator, I interpolation, and A averaging
% the reduced form is
% B' * Cx * ( (I*B*R + yI).*(A*B*R + yA) ) =
% which is split into four components:
% B' * Cx * ( (I*B).*(A*B) ) * kron(R,R) +  [quadratic in R]
% B' * Cx * ( (I*B).*(yA) + yI.*(A*B) ) * R +  [linear in R]
% B' * Cx * ( yA.*yI ) [constant]

% if we include Vbc:
% Cx * ( (I*(B*R+Vbc) + yI).*(A*(B*R+Vbc) + yA) ) = 
% Cx * ( (I*B*R + I*Vbc + yI).*(A*B*R + A*Vbc + yA) ) = 
% Cx * ( (I*B*R).*(A*B*R) +
%        (I*B*R).*(A*Vbc + yA) + (I*Vbc + yI).*(A*B*R) +
%        (I*Vbc+yI).*(A*Vbc+yA) ) 

% if we evaluate convection with Vbc as argument, we get the last term:
% Cx * ( (I*Vbc + yI).*(A*Vbc + yA) )

% note: first argument of convection is V, second is C
% conv(Vbc,Vbc)
[convu, convv] = convection(Vbc,Vbc,0,options,0);
conv_bc        = B'*(Diag.*[convu;convv]);

% the following calls are needed for constructing conv_linear
% conv(0,0)
[convu, convv] = convection(Z,Z,0,options,0);
conv_bc0       = B'*(Diag.*[convu;convv]);
% conv(Vbc,0)
[convu, convv] = convection(Vbc,Z,0,options,0);
conv_bc1       = B'*(Diag.*[convu;convv]);
% conv(0,Vbc)
[convu, convv] = convection(Z,Vbc,0,options,0);
conv_bc2       = B'*(Diag.*[convu;convv]);

% second, compute the linear terms
conv_linear1 = zeros(M,M);
conv_linear2 = zeros(M,M);
conv_linear3 = zeros(M,M);
conv_linear4 = zeros(M,M);

for i=1:M
    % we need to evaluate
    % (1): (I*B*R).*(A*Vbc + yA) 
    % and
    % (2): (I*Vbc + yI).*(A*B*R)

    % first evaluate conv(Vbc,B):
    % Cx * ( (I*B + yI).*(A*Vbc + yA) )
    [convu1, convv1]   = convection(Vbc,B(:,i),0,options,0);
    % and conv(B,Vbc):
    % Cx * ( (I*Vbc + yI).*(A*B + yA) )
    [convu2, convv2]   = convection(B(:,i),Vbc,0,options,0);
    
    % linear term (1) should be
    % Cx * ( (I*B).*(A*Vbc + yA) ) 
    % which can be computed from
    % Cx * ( (I*B + yI).*(A*Vbc + yA) ) - Cx * ( yI.*(A*Vbc + yA) )
    % conv(Vbc,B) - conv(Vbc,0)
    conv_linear1(:,i)  = B'*(Diag.*[convu1;convv1]) - conv_bc1;
    
    % in a similar fashion the second linear term:
    % Cx * ( (I*Vbc + yI).*(A*B) ) 
    % = Cx * ( (I*Vbc + yI).*(A*B + yA) ) - Cx * ( (I*Vbc + yI).*yA )
    % conv(B,Vbc) - conv(0,Vbc)
    conv_linear2(:,i)  = B'*(Diag.*[convu2;convv2]) - conv_bc2;
    
    % get conv(0,B) to be used for quadratic terms
    [convu3, convv3]   = convection(Z,B(:,i),0,options,0);
    conv_linear3(:,i)  = B'*(Diag.*[convu3;convv3]);
    % get conv(B,0) for quadratic terms
    [convu4, convv4]   = convection(B(:,i),Z,0,options,0);
    conv_linear4(:,i)  = B'*(Diag.*[convu4;convv4]);
end

% quadratic terms
% this can be done in different ways; e.g. a third order tensor of size (M,M,M)
% or a matricized tensor of size M*(M^2) or (M^2)*M
% here we choose for M*(M^2), like in Kramer & Willcox (AIAA 2019)

% we want to get 
% Cx ( (I*B*R).*(A*B*R) ) = 
% Cx * ( (I*B*R + yI).*(A*B*R + yA) ) - Cx * ( (I*B*R+yI) .* yA  + yI.*(A*B*R +
% yA) ) + Cx*(yI.*yA)
% = conv(B,B) - conv(0,B) - conv(B,0) + conv(0,0)
conv_quad   = zeros(M,M*M);
for i=1:M
    for j=1:M
        [convu, convv] = convection(B(:,i),B(:,j),0,options,0);
        % convert third order tensor to a matrix via the following indexing:
        k = sub2ind([M,M],j,i);
        conv_quad(:,k) = B'*(Diag.*[convu;convv]) - conv_linear3(:,i) ...
                          - conv_linear4(:,j) + conv_bc0;        
    end
end

options.rom.Conv_quad   = conv_quad;
options.rom.Conv_linear = conv_linear1 + conv_linear2;
options.rom.yConv       = conv_bc;
