function options = operator_rom(Vbc,options)

B  = options.rom.B;

if (options.rom.weighted_norm == 0)
    Diag = options.grid.Om_inv;
elseif (options.rom.weighted_norm == 1)
    Diag = ones(options.grid.Nu+options.grid.Nv,1);
end

% we construct the ROM operators in the following 'non-intrusive' way, i.e.
% we hardly require any knowledge of the exact FOM discretization. we
% simply make calls to diffusion and convection with the modes phi as input
% arguments

% for diffusion:
% the ROM discretization will read Diff*R + yDiff, where Diff is the reduced
% diffusion operator
% we get this by multiplying B'*Diffusion(B*R), where a call to Diffusion(V) returns
% D*V + yD = D*(B*R) + yD
M    = options.rom.M;
Diff = zeros(M,M);

% first evaluate diffusion with zero velocity (or Vbc field) to only get boundary
% condition contribution
[d2u_bc,d2v_bc] = diffusion(Vbc,0,options,0);
yDiff           = B'*(Diag.*[d2u_bc;d2v_bc]);

% now evaluate each column of the reduced matrix as B'*D*B, where the
% boundary conditions are subtracted in order to get only the terms
% that involve the matrix
for i=1:M
    [d2u,d2v] = diffusion(B(:,i),0,options,0);
    Diff(:,i) = B'*(Diag.*[d2u;d2v]) - yDiff;
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

% first compute the constant term, by evaluating convection with zero velocities
% (or Vbc field) to only get
[convu, convv] = convection(Vbc,Vbc,0,options,0);
conv_bc        = B'*(Diag.*[convu;convv]);

% second, compute the linear terms
conv_linear1 = zeros(M,M);
conv_linear2 = zeros(M,M);
for i=1:M
    [convu1, convv1]   = convection(B(:,i),Vbc,0,options,0);
    [convu2, convv2]   = convection(Vbc,B(:,i),0,options,0);
    
    conv_linear1(:,i)  = B'*(Diag.*[convu1;convv1]) - conv_bc;
    conv_linear2(:,i)  = B'*(Diag.*[convu2;convv2]) - conv_bc;
end

% quadratic terms
% this can be done in different ways; e.g. a third order tensor of size (M,M,M)
% or a matricized tensor of size M*(M^2) or (M^2)*M
% here we choose for M*(M^2), like in Kramer & Willcox (AIAA 2019)
conv_quad   = zeros(M,M*M);
for i=1:M
    for j=1:M
        [convu, convv] = convection(B(:,i),B(:,j),0,options,0);
        % convert third order tensor to a matrix via the following indexing:
        k = sub2ind([M,M],j,i);
        conv_quad(:,k) = B'*(Diag.*[convu;convv]) - conv_linear1(:,i) ...
                          - conv_linear2(:,j) - conv_bc;        
    end
end

options.rom.Conv_quad   = conv_quad;
options.rom.Conv_linear = conv_linear1 + conv_linear2;
options.rom.yConv       = conv_bc;
