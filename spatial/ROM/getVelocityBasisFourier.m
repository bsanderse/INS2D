function [B,Bp,M] = getVelocityBasisFourier(options)
        
%% Fourier
Nu  = options.grid.Nu;
Nv  = options.grid.Nv;

% NV  = options.grid.NV;

% assume homogeneous boundary conditions, so no Vbc field and no yM
% Vbc    = zeros(NV,1);
% rom_yM = 0;


M = options.rom.M;

% construct velocity basis B from discrete Fourier transform
% note that we have V = B*R, which we split into two components
% u = B_u*R_u, v = B_v*Rv
Phi = getFourierBasis(options); %,M);
% scale with sqrt(Omega) in order to get Phi'*Om*Phi = I
Phi = 1/sqrt(options.grid.Om(1)) * Phi;
B   = [Phi spalloc(Nu,M,0); ...
       spalloc(Nv,M,0) Phi];

% M can be different from the input given by options.rom.M
M   = size(B,2);
% options.rom.M = M;

% we also need to make a pressure basis Bp because B is not
% divergence-free
%         Bp = B;
Bp = Phi;

% set pressure_recovery to 2 in order that reduced Poisson operator
% is constructed
div_basis = max(abs(options.discretization.M*B),[],1); %
% max over all columns:
maxdiv_basis = max(div_basis);
if (maxdiv_basis > 1e-14)
    disp(['Fourier basis is not divergence free: ' num2str(maxdiv_basis)]);
    disp('Adding basis for pressure');
%     div_free = 0;
% else     
%     div_free = 1;
end

%        
%         % make the basis divergence free
%         Om_inv = options.grid.Om_inv;
%         G  = options.discretization.G;
%         for i=1:M
%            
%             % make velocity field divergence free            
%         
%             f  = options.discretization.M*B(:,i) + options.discretization.yM;
%             dp = pressure_poisson(f,t,options);
%             
%             % look at the resulting modes, before and after projection
% %             figure(1)
% %             surf(reshape(B(options.grid.indu,i),Nx,Ny)')
%             
%             B(:,i)  = B(:,i) - Om_inv.*(G*dp);
%             % turns out that they effectively become 0!
% %             figure(2)
% %             surf(reshape(B(options.grid.indu,i),Nx,Ny)')
%         end
%         
%         % make orthogonal
%         [Bnew,R] = qr(B);
% %         B = Bnew/sqrt(options.grid.Om(1));
%         
%         [B,R] = mwgson(B,options.grid.Om);

%         div_basis = max(abs(options.discretization.M*B),[],1); %
%         % max over all columns:
%         maxdiv_basis = max(div_basis);
%         if (maxdiv_basis > 1e-14)
%             disp(['Fourier basis is not divergence free: ' num2str(maxdiv_basis)]);
%             disp('Adding basis for pressure');
%             options.rom.div_free = 0;
%         else     
%             options.rom.div_free = 1;
%         end

% we do not expect orthogonality but we do expect a diagonal structure
disp('orthogonality of reduced basis B wrt omega:');
Om     = options.grid.Om;
Om_mat = spdiags(Om,0,Nu+Nv,Nu+Nv);
B_orth = max(max(abs(B'*Om_mat*B-speye(size(B,2)))))
if (B_orth>1e-13)
    warning('Reduced basis B not orthogonal wrt Omega');
end

keyboard