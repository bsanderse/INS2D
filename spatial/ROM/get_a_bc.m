function a_bc = get_a_bc(t,options)
if options.rom.BC_DEIM == 0
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

else
    %% DEIM botch
%     y_bc = get_bc_vector_yBC(t,options);
%     PTf = y_bc(options.rom.bc_deim_inds);
%     a_bc = doDEIM(options.rom.bc_deim_PTUinv,PTf);

    %% testing
%     if options.rom.BC_DEIMdim == options.rom.Mbc
%         norm(a_bc - get_a_bc_precomp(t,options))
%     end

% efficient
    PTybc = ybc_entries(options.rom.bc_deim_inds, ...
        options.rom.ybc_coords,t,options);
    a_bc = doDEIM(options.rom.bc_deim_PTUinv,PTybc);

% inefficient
%     PTybc2 = ybc_entries2(options.rom.bc_deim_inds, ...
%         options.rom.ybc_coords,t,options);


%     a_bc_effic = doDEIM(options.rom.bc_deim_PTUinv,PTybc);
% %     norm(a_bc - a_bc_effic)
%     norm(PTybc - PTybc2)
end