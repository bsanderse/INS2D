function [conv_bc,conv_linear,conv_quad] = operator_rom_convection(P,options)
% precompute convective operators
% projection with generic matrix P, size M1 x NV
% for momentum equation, P is B' or B'*Om_inv
% for Poisson equation, P is Bp'*M*Om_inv

NV = options.grid.Nu+options.grid.Nv;
Z  = zeros(NV,1);

B   = options.rom.B;
M   = options.rom.M;
Vbc = options.rom.Vbc;

% P can project to M or to Mp modes 
M1  = size(P,1);   

%% convection
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
conv_bc        = P*[convu;convv];

% the following calls are needed for constructing conv_linear
% conv(0,0)
[convu, convv] = convection(Z,Z,0,options,0);
conv_bc0       = P*[convu;convv];
% conv(Vbc,0)
[convu, convv] = convection(Vbc,Z,0,options,0);
conv_bc1       = P*[convu;convv];
% conv(0,Vbc)
[convu, convv] = convection(Z,Vbc,0,options,0);
conv_bc2       = P*[convu;convv];

% second, compute the linear terms
conv_linear  = zeros(M1,M);
conv_linear3 = zeros(M1,M);
conv_linear4 = zeros(M1,M);

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
%     conv_linear1(:,i)  = B'*(Diag.*[convu1;convv1]) - conv_bc1;
    
    % in a similar fashion the second linear term:
    % Cx * ( (I*Vbc + yI).*(A*B) ) 
    % = Cx * ( (I*Vbc + yI).*(A*B + yA) ) - Cx * ( (I*Vbc + yI).*yA )
    % conv(B,Vbc) - conv(0,Vbc)
%     conv_linear2(:,i)  = B'*(Diag.*[convu2;convv2]) - conv_bc2;
    
    % directly create conv_linear, which is conv_linear1+conv_linear2
    conv_linear(:,i)   = P*[convu1;convv1] - conv_bc1 + ...
                         P*[convu2;convv2] - conv_bc2;
    
    % get conv(0,B) to be used for quadratic terms
    [convu3, convv3]   = convection(Z,B(:,i),0,options,0);
    conv_linear3(:,i)  = P*[convu3;convv3];
    % get conv(B,0) for quadratic terms
    [convu4, convv4]   = convection(B(:,i),Z,0,options,0);
    conv_linear4(:,i)  = P*[convu4;convv4];
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
conv_quad   = zeros(M1,M*M);
for i=1:M
    for j=1:M
        [convu, convv] = convection(B(:,i),B(:,j),0,options,0);
        % convert third order tensor to a matrix via the following indexing:
        k = sub2ind([M,M],j,i); % first loop over j,  then over i
        conv_quad(:,k) = P*[convu;convv] - conv_linear3(:,i) ...
                          - conv_linear4(:,j) + conv_bc0;        
    end
end

% options.rom.Conv_quad   = conv_quad;
% options.rom.Conv_linear = conv_linear1 + conv_linear2;
% options.rom.yConv       = conv_bc;