
if (xfoil==1)
    %% load from Xfoil data

    % on the mac: (on linux this is probably simply 'xfoil')
    xfoil_path   = '/Users/sanderse/Software/Xfoil/bin/xfoil';

    ext          = '.dat';
    
    xfoil_folder     = 'force/airfoil_data/';
    addpath(xfoil_folder);
    xfoil_file       = [xfoil_folder 'naca_' num2str(airfoil_type) '_' num2str(aoa) '_' num2str(Re_c)];
    xfoil_file_input = [xfoil_file '_input' ext];
    xfoil_file_cp    = [xfoil_file '_cp' ext];
    xfoil_file_cf    = [xfoil_file '_cf' ext]; 
    xfoil_file_coord = [xfoil_file '_coord' ext];
    xfoil_file_camb  = [xfoil_file '_camb' ext];

    
    % write a file with xfoil commands
        fid = fopen(xfoil_file_input, 'w');
        if fid == -1
            error('Cannot open file for writing.');
        end

        nl = sprintf('\n');
        % naca
        fwrite(fid,['naca ' airfoil_type nl]);
        % close trailing edge with blending 0.2
        fwrite(fid,['gdes ' nl]);
        fwrite(fid,['gset ' nl]);
        fwrite(fid,['tgap 0 0.2' nl]);
%         fwrite(fid,['adeg 10' nl]);
        fwrite(fid,['exec ' nl]);  
        % write camber line line
        fwrite(fid,['camb ' nl]);
        fwrite(fid,['wrtc ' nl]);
        fwrite(fid,[xfoil_file_camb nl]);
        fwrite(fid,nl);            
        fwrite(fid,nl);
        fwrite(fid,nl);
        % save coordinates
        fwrite(fid,['psav ' xfoil_file_coord nl]);
        fwrite(fid,nl);
        % operating menu
        fwrite(fid,['oper' nl]);
        if (Re_c~=0)
            % set Reynolds
            fwrite(fid,['visc' nl]);
            fwrite(fid,[num2str(Re_c) nl]);
        end
        % set angle of attack
        fwrite(fid,['alfa ' num2str(aoa) nl]);
        % write Cp
        fwrite(fid,['cpwr ' xfoil_file_cp nl]);
        if (Re_c~=0)
            % write boundary layer data
            fwrite(fid,['dump ' xfoil_file_cf nl]);
        end
        fwrite(fid,nl);
        fwrite(fid,'quit');

%         system([xfoil_path ' <' num2str(xfoil_file_input) '> /dev/null']);

        % read xfoil data in matlab
        [hdr1, airfoil_Cp] = hdrload(xfoil_file_cp);
        [hdr2, airfoil_Cf] = hdrload(xfoil_file_cf);
        airfoil_coord = load(xfoil_file_coord);
%         [hdr3, airfoil_camb] = hdrload(xfoil_file_camb);


        % coordinates are panel endpoints
        x_k  = squeeze(airfoil_coord(:,1));
        y_k  = squeeze(airfoil_coord(:,2));
%         x_cl = squeeze(airfoil_camb(:,1));
%         y_cl = squeeze(airfoil_camb(:,2));
        
                
        % number of endpoints is nk, number of panels nk-1:
        % the begin and end point are equal
        nk   = length(x_k);
        if (rem(nk,2)==0) % nk even
            nk2  = nk/2;
        else
            disp('number of points should be even');
        end        
        
        % Cp at panel endpoints
        Cp_ep   = squeeze(airfoil_Cp(:,2));
        xCp_ep  = squeeze(airfoil_Cp(:,1));
        if (Re~=0)
            Cf_ep   = squeeze(airfoil_Cf(:,7));
            sCf_ep  = squeeze(airfoil_Cf(:,1));
            xCf_ep  = squeeze(airfoil_Cf(:,2));
            
            % cut off the 'wake' part
            Cf_ep   = Cf_ep(1:nk);
            sCf_ep  = sCf_ep(1:nk);
        else
            Cf_ep   = zeros(size(Cp_ep));
        end

        if (max(abs(xCp_ep-x_k))>1e-3)
            warning('difference between xk and xCp coordinates too large');
        end
       
        % panel midpoints ('collocation points')
        % nk-1 points (closed trailing edge)
        [xcol,ycol,Sk,vect,vecn] = panel_properties(x_k,y_k,1,1);
%         plotnormal(x_k,y_k,vect,vecn);
        
        % define leading edge as point with minimum x_k; this point belongs
        % to both upper and lower side 
        [val ile] = min(x_k);
        % x_k(ile+1) can also be a minimum (e.g. symmetric airfoil)
        
        % get Cp at panel midpoints (collocation points)
        % assume nk2 points upper side, nk2 points lower side
        Cp_u = interp1(x_k(1:ile),Cp_ep(1:ile),xcol(1:ile-1),'pchip','extrap');
        Cp_l = interp1(x_k(ile:end),Cp_ep(ile:end),xcol(ile:end),'pchip','extrap');
        Cp   = [Cp_u;Cp_l];
        
        % get Cf at panel midpoints (collocation points)
        Cf_u = interp1(x_k(1:ile),Cf_ep(1:ile),xcol(1:ile-1),'pchip','extrap'); % spline gives Gibbs phenomena near discontinuities
        Cf_l = interp1(x_k(ile:end),Cf_ep(ile:end),xcol(ile:end),'pchip','extrap'); 
        Cf   = [Cf_u;Cf_l];        
        
   
else
    %% own panel method
    nk        = 200; 
    nk2       = floor(nk/2);
    i         = 1:nk2+1;
    
    % clockwise numbering from LE to TE and back
    x_kb      = sin(0.5*pi*(i-1)*(1/nk2)).^2;   
    x_ko      = sin(0.5*pi*(1-(i-1)*(1/nk2))).^2;  

%     if (~strcmp(airfoil_type,'0012'))
%         error('not implemented');
%     end
    if (length(airfoil_type)>4)
        error('not implemented');
    end
    naca_m = str2double(airfoil_type(1))/100;
    naca_p = str2double(airfoil_type(2))/10;
    naca_t = str2double(airfoil_type(3:4))/100;
    
    % naca thickness without camber:
    y_kb_t      = naca_t/0.2*(0.29690*x_kb.^0.5 - 0.12600*x_kb - 0.35160*x_kb.^2 + 0.28430*x_kb.^3 -0.10150*x_kb.^4);
%     y_ko_t      = -naca_t/0.2*(0.29690*x_ko.^0.5 - 0.12600*x_ko - 0.35160*x_ko.^2 + 0.28430*x_ko.^3 -0.10150*x_ko.^4);
    % naca camberline
    if (naca_m~=0 && naca_p~=0)
        y_c      = naca_m/(naca_p^2)*x_kb.*(2*naca_p-x_kb).*(x_kb<=naca_p) + ...
                   naca_m/(1-naca_p)^2*(1-x_kb).*(1+x_kb-2*naca_p).*(x_kb>naca_p);
             
        dycdx =  2*naca_m/(naca_p^2) * (naca_p-x_kb).*(x_kb<=naca_p) + ...
                 2*naca_m/((1-naca_p)^2) * (naca_p-x_kb).*(x_kb>naca_p);
             
        theta = atan(dycdx);
        x_kb  = x_kb - y_kb_t.*sin(theta);
        y_kb  = y_c + y_kb_t.*cos(theta);
        x_ko  = x_kb + y_kb_t.*sin(theta);
        y_ko  = y_c - y_kb_t.*cos(theta);
        x_ko  = fliplr(x_ko);
        y_ko  = fliplr(y_ko);
    end
    
    
    % cylinder
    % y_kb      =  sqrt(0.25-(x_kb-0.5).^2);
    % y_ko      = -sqrt(0.25-(x_ko-0.5).^2);

    x_k       = [x_kb(1:end-1) x_ko(1:end)];
    y_k       = [y_kb(1:end-1) y_ko(1:end)];

    % find Cp with panel method
    aoa       = 0;
    [Cp,Vs,xcol,ycol,Sk,vecn,vect] = Panel(x_k,y_k,aoa);    % call panel method to calculate Cp distribution
        
    % camber line
    % assume that upper and lower distribution have same x coordinates
    x_cl      = x_kb(1:end);
    y_cl      = 0.5*(y_kb(1:end)+y_ko(end:-1:1));
    % Cf is not known
    Cf        = zeros(size(Cp));
end


CpSx  = Cp.*Sk.*vecn(:,1);
CfSx  = Cf.*Sk.*vect(:,1);
CpSy  = Cp.*Sk.*vecn(:,2);
CfSy  = Cf.*Sk.*vect(:,2);

% a negative Cp is in the direction of the unit vector
% assuming the chord is aligned with the x-axis, as is the case normally
% with xfoil computations
Cx    = -CpSx+CfSx;
Cy    = -CpSy+CfSy;

% total forces
CX    = sum(Cx);
CY    = sum(Cy);

% local force in lift and drag directions
% Cl    = -Cx*sind(aoa) + Cy*cosd(aoa);
% Cd    = Cx*cosd(aoa) + Cy*sind(aoa);

% as a check with Xfoil
% Xfoil calculates Cl from contour integration, but Cd is computed from
% wake analysis (Squire-Young). obtaining CD from contour integration gives
% large errors, especially due to the contribution of Cp to the x
% direction.
% note furthermore that errors have been introduced in the interpolation from endpoints
% to midpoints
Cl    = -CX*sind(aoa) + CY*cosd(aoa)
Cd    = CX*cosd(aoa) + CY*sind(aoa)
