function a_bc = get_a_bc(t,options)

% if abs(options.rom.time_vec - t) > 1e-12
%     error("a_bc dictionary does not work properly")
% end
% 
% a_bc = options.rom.a_bc_vec;

%% botch

a_bc = get_a_bc_precomp(t,options); %

%     phi_bc = options.rom.phi_bc;
%     
%     a_bc = phi_bc'*get_bc_vector_yBC(t,options);