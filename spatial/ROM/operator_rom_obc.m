function [obc_hom,obc_inhom,obc_hom_inhom,obc_inhom_hom,obc_hom_inhom2] ...
    = operator_rom_obc(P,options)

M1  = size(P,1);   % P =/= B' possible

B   = options.rom.B;
M   = options.rom.M;
phi_inhom = options.rom.phi_inhom;
M_inhom = size(phi_inhom,2);

Conv_diag = options.grid.C;

obc_hom = zeros(M1,M*M);
for i=1:M
    for j=1:M
        obc = (Conv_diag*B(:,i)).*(B(:,j));
        % convert third order tensor to a matrix via the following indexing:
        k = sub2ind([M,M],j,i); % first loop over j,  then over i
        obc_hom(:,k) = P*obc;        
    end
end

obc_inhom = zeros(M1,M_inhom*M_inhom);
for i=1:M_inhom
    for j=1:M_inhom
        obc = (Conv_diag*phi_inhom(:,i)).*(phi_inhom(:,j));
        % convert third order tensor to a matrix via the following indexing:
        k = sub2ind([M_inhom,M_inhom],j,i); % first loop over j,  then over i
        obc_inhom(:,k) = P*obc;        
    end
end

%% use obc_hom_inhom2 to instead
% obc_hom_inhom = zeros(M1,M*M_inhom);
% for i=1:M
%     for j=1:M_inhom
%         obc = (Conv_diag*B(:,i)).*(phi_inhom(:,j));
%         % convert third order tensor to a matrix via the following indexing:
%         k = sub2ind([M_inhom,M],j,i); % first loop over j,  then over i
%         obc_hom_inhom(:,k) = P*obc;        
%     end
% end
% 
% obc_inhom_hom = zeros(M1,M_inhom*M);
% for i=1:M_inhom
%     for j=1:M
%         obc = (Conv_diag*phi_inhom(:,i)).*(B(:,j));
%         % convert third order tensor to a matrix via the following indexing:
%         k = sub2ind([M,M_inhom],j,i); % first loop over j,  then over i
%         obc_inhom_hom(:,k) = P*obc;        
%     end
% end
%%
obc_hom_inhom2 = zeros(M1,M*M_inhom);
for i=1:M
    for j=1:M_inhom
        obc = (Conv_diag*B(:,i)).*(phi_inhom(:,j)) ...
            + (Conv_diag*phi_inhom(:,j)).*(B(:,i));
        % convert third order tensor to a matrix via the following indexing:
        k = sub2ind([M_inhom,M],j,i); % first loop over j,  then over i
        obc_hom_inhom2(:,k) = P*obc;        
    end
end