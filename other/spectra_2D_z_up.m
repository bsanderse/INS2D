% Does 1D FFT along flow directed grid lines at specific height and then
% averaged in spanwise direction.
function [] = spectra_2D_z_up(dir,plot_dir,cases,variable_name,data_type,binary_type,zipped,save_on,print_on,nx,ny,nz,x_min,y_min,z_min,x_max,y_max,z_max,heights);



%% Figure vertical indices corresponding to given heights
nh = length(heights);

fidh = fopen([dir,'/averaging/0/hlevelsCell'],'r');
z = fscanf(fidh,'%g',inf);
fclose(fidh);

for i = 1:nh
    [heights_act(i) heights_index(i)] = min(abs(z-heights(i)));
    heights_act(i) = z(heights_index(i));
end


%% Read in the data
nt = length(cases);

if (nx <= ny)
    ne = floor(nx/2)+1;
else
    ne = floor(ny/2)+1;
end
e1plot = zeros(nh,ne);
e2plot = zeros(nh,ne);
e3plot = zeros(nh,ne);
e12plot = zeros(nh,ne);
e123plot = zeros(nh,ne);

for n = 1:nt
    disp(['Processing case: ',num2str(cases(n)),'...']);
    for h = 1:nh
        tic
        disp(['   Desired Height: ',num2str(heights(h)),' m']);
        disp(['   Actual Height:  ',num2str(heights_act(h)),' m']);
        
        % Open the file
        if (zipped == 1)
            gunzip([dir,'/',num2str(cases(n)),'/',variable_name,'.gz']);
        end
        fid = fopen([dir,'/',num2str(cases(n)),'/',variable_name],'r');
        
        % Check to see if binary (binary_flag = 1) or ascii (binary_flag = 0)
        flag = 0;
        while (flag == 0)
            data = fgetl(fid);
            k = findstr(data,'format');
            if (isempty(k) == 0)
                flag = 1;
                k = findstr(data,'ascii');
                if (isempty(k) == 0)
                    binary_flag = 0;
                    disp('   Data is ascii');
                end
                k = findstr(data,'binary');
                if (isempty(k) == 0)
                    binary_flag = 1;
                    disp('   Data is binary');
                end
            end
        end
        
        % Find the end of the header
        flag = 0;
        while (flag == 0)
            data = fgetl(fid);
            k = findstr(data,'internalField');
            if (isempty(k) == 0)
                flag = 1;
            end
        end
        
        % Find the number of cells
        data = fgetl(fid);
        cells = str2num(data);
        
        % Read in the data
        % ascii
        if (binary_flag == 0)
            % read up to specified position in file
            data = fgetl(fid);
            if     (data_type == 1)
                fscanf(fid,'%g\n',[1 nx*ny*(heights_index(h)-1)])';
            elseif (data_type == 2)
                fscanf(fid,'(%g %g %g)\n',[3 nx*ny*(heights_index(h)-1)])';
            elseif (data_type == 3)
                fscanf(fid,'(%g %g %g %g %g %g)\n',[6 nx*ny*(heights_index(h)-1)])';
            elseif (data_type == 4)
                fscanf(fid,'(%g %g %g %g %g %g %g %g %g)\n',[9 nx*ny*(heights_index(h)-1)])';
            end
            
            % read data at specified height
            if     (data_type == 1)
                var = fscanf(fid,'%g\n',[1 nx*ny])';
            elseif (data_type == 2)
                var = fscanf(fid,'(%g %g %g)\n',[3 nx*ny])';
            elseif (data_type == 3)
                var = fscanf(fid,'(%g %g %g %g %g %g)\n',[6 nx*ny])';
            elseif (data_type == 4)
                var = fscanf(fid,'(%g %g %g %g %g %g %g %g %g)\n',[9 nx*ny])';
            end
            
            % binary
        elseif (binary_flag == 1)
            % read up to specified position in file
            data = fread(fid,1,'char');
            if     (data_type == 1)
                fseek(fid,1*nx*ny*(heights_index(h)-1)*8,'cof');
            elseif (data_type == 2)
                fseek(fid,3*nx*ny*(heights_index(h)-1)*8,'cof');
            elseif (data_type == 3)
                fseek(fid,6*nx*ny*(heights_index(h)-1)*8,'cof');
            elseif (data_type == 4)
                fseek(fid,9*nx*ny*(heights_index(h)-1)*8,'cof');
            end
            
            % read data at specified height
            if     (data_type == 1)
                var = fread(fid,[1 nx*ny],binary_type)';
            elseif (data_type == 2)
                var = fread(fid,[3 nx*ny],binary_type)';
            elseif (data_type == 3)
                var = fread(fid,[6 nx*ny],binary_type)';
            elseif (data_type == 4)
                var = fread(fid,[9 nx*ny],binary_type)';
            end
            
        end
        
        % Close the file
        fclose(fid);
        
        
        %% OpenFOAM data is unstructured, so structure the data
        U = zeros(nx,ny);
        V = zeros(nx,ny);
        W = zeros(nx,ny);
        
        for j = 1:ny
            U(1:nx,j) = var(((j-1)*nx)+1:(j*nx),1);
            V(1:nx,j) = var(((j-1)*nx)+1:(j*nx),2);
            W(1:nx,j) = var(((j-1)*nx)+1:(j*nx),3);
        end
        
        
        %% Take averages and find fluctuations
        Ubar = mean(mean(U));
        Vbar = mean(mean(V));
        Wbar = mean(mean(W));
        Uprime = U-Ubar;
        Vprime = V-Vbar;
        Wprime = W-Wbar;
        
        k_p = 0.5*(sum(sum(Uprime.*Uprime)) + sum(sum(Vprime.*Vprime)) + sum(sum(Wprime.*Wprime)))/(nx*ny);
        
        
        %% Perform FFT
        delta_x = (x_max-x_min)/nx;
        delta_y = (y_max-y_min)/ny;
        x = x_min+(delta_x/2.0):delta_x:x_max-(delta_x/2.0);
        y = y_min+(delta_y/2.0):delta_y:y_max-(delta_y/2.0);
        
        %  set x and y wavenumbers in reverse wrap-around format (don't scale yet,
        %  just integers)
        k1 = zeros(nx,1);
        k2 = zeros(ny,1);
        k1(1:floor(nx/2)+1) = (0:1:floor(nx/2));
        k1(floor(nx/2)+2:nx) = (-floor(nx/2)+1:1:-1);
        k2(1:floor(ny/2)+1) = (0:1:floor(ny/2));
        k2(floor(ny/2)+2:ny) = (-floor(ny/2)+1:1:-1);
        
        %  make a matrix of the wavenumbers
        [K1,K2] = meshgrid(k1,k2);
        
        %  make a matrix of the magnitude of the wavevector made up of the kx and
        %  ky components
        K12 = sqrt(K1.^2 + K2.^2)';
        
        %  make a wavenumber magnitude list corresponding to the half range of the
        %  smaller of the two dimensions, x or y.
        if (nx <= ny)
            k2d = k1(1:floor(nx/2)+1);
        else
            k2d = k2(1:floor(ny/2)+1);
        end
        
        %  initialize the spectra
        e1 = zeros(size(k2d))';
        e2 = zeros(size(k2d))';
        e3 = zeros(size(k2d))';
        e12 = zeros(size(k2d))';
        e123 = zeros(size(k2d))';
        
        %  take the 2d fft, and the sum the energy within discrete ranges of the
        %  wavevector magnitude.
        Uhat = fft2(Uprime);
        for k=1:length(k2d)
            e1(k) = 0.0;
            count = 0;
            for i = 1:nx
                for j = 1:ny
                    if ( (K12(i,j) < k2d(k)+0.5) && (K12(i,j) > k2d(k)-0.5) )
                        count = count+1;
                        e1(k) = e1(k)+abs(Uhat(i,j))^2;
                    end
                end
            end
            %eu(m,k) = eu(m,k)/count;
        end
        
        Vhat = fft2(Vprime);
        for k=1:length(k2d)
            e2(k) = 0.0;
            count = 0;
            for i = 1:nx
                for j = 1:ny
                    if ( (K12(i,j) < k2d(k)+0.5) && (K12(i,j) > k2d(k)-0.5) )
                        count = count+1;
                        e2(k) = e2(k)+abs(Vhat(i,j))^2;
                    end
                end
            end
            %ev(m,k) = ev(m,k)/count;
        end
        
        What = fft2(Wprime);
        for k=1:length(k2d)
            e3(k) = 0.0;
            count = 0;
            for i = 1:nx
                for j = 1:ny
                    if ( (K12(i,j) < k2d(k)+0.5) && (K12(i,j) > k2d(k)-0.5) )
                        count = count+1;
                        e3(k) = e3(k)+abs(What(i,j))^2;
                    end
                end
            end
            %ew(m,k) = ew(m,k)/count;
        end      
        
        UVhat = abs(Uhat).^2+abs(Vhat).^2;
        for k=1:length(k2d)
            e12(k) = 0.0;
            count = 0;
            for i = 1:nx
                for j = 1:ny
                    if ( (K12(i,j) < k2d(k)+0.5) && (K12(i,j) > k2d(k)-0.5) )
                        count = count+1;
                        e12(k) = e12(k)+UVhat(i,j);
                    end
                end
            end
            %euvw(m,k) = euvw(m,k)/count;
        end
        
        UVWhat = abs(Uhat).^2+abs(Vhat).^2+abs(What).^2;
        for k=1:length(k2d)
            e123(k) = 0.0;
            count = 0;
            for i = 1:nx
                for j = 1:ny
                    if ( (K12(i,j) < k2d(k)+0.5) && (K12(i,j) > k2d(k)-0.5) )
                        count = count+1;
                        e123(k) = e123(k)+UVWhat(i,j);
                    end
                end
            end
            %euvw(m,k) = euvw(m,k)/count;
        end
        
        e1plot(h,:)   = e1plot(h,:) + (e1./((nx*ny)^2));
        e2plot(h,:)   = e2plot(h,:) + (e2./((nx*ny)^2));
        e3plot(h,:)   = e3plot(h,:) + (e3./((nx*ny)^2));
        e12plot(h,:)  = e12plot(h,:)+(e12./((nx*ny)^2));
        e123plot(h,:) = e123plot(h,:)+(e123./((nx*ny)^2));
        
        k_s2 = 0.5*(sum(sum(UVWhat))./((nx*ny)^2));
        k_diff = 100*(k_s2-k_p)/k_p;
        elapsed_time = toc;
        disp(['   Energy (physical space: ',num2str(k_p),' (m/s)^2']);
        disp(['   Energy (spectral space: ',num2str(k_s2),' (m/s)^2']);
        disp(['   Difference:             ',num2str(k_diff),' %']);
        disp(['   Elapsed time:           ',num2str(elapsed_time),' s']);
        disp(' ');
    end
end
% e1plot = e1plot(:,1:floor(nx/2)+1)./(ny*nt);
% e2plot = e2plot(:,1:floor(nx/2)+1)./(ny*nt);
% e3plot = e3plot(:,1:floor(nx/2)+1)./(ny*nt);
% e123plot = e123plot(:,1:floor(nx/2)+1)./(ny*nt);


%% Scale the wavevector magnitude by the coarsest grid spacing
if (nx <= ny)
    k2d = k2d*(2.0*pi/(nx*delta_x));
    kappa_max = (2.0*pi/(2.0*delta_x));
    kappa_max_23 = (2.0/3.0)*kappa_max;
else
    k2d = k2d*(2.0*pi/(ny*delta_y));
    kappa_max = (2.0*pi/(2.0*delta_y));
    kappa_max_23 = (2.0/3.0)*kappa_max;
end


%% Save data
if (save_on == 1)
    
    if (exist([dir,'/matlab_data']) == 0)
        mkdir([dir,'/'],'matlab_data')
    end
    if (exist([dir,'/spectra_data']) == 0)
        mkdir([dir,'/'],'spectra_data')
    end
    
    save([dir,'/matlab_data/spectra_2d.mat'],'heights','heights_act','k2d','e1plot','e2plot','e3plot','e12plot','e123plot');
    
    for h = 1:nh
        fid1   = fopen([dir,'/spectra_data/e1x_',num2str(heights(h)),'.dat'],'w');
        fid2   = fopen([dir,'/spectra_data/e2x_',num2str(heights(h)),'.dat'],'w');
        fid3   = fopen([dir,'/spectra_data/e3x_',num2str(heights(h)),'.dat'],'w');
        fid12 = fopen([dir,'/spectra_data/e12x_',num2str(heights(h)),'.dat'],'w');
        fid123 = fopen([dir,'/spectra_data/e123x_',num2str(heights(h)),'.dat'],'w');
        nk = length(k2d);
        for i = 1:nk
            fprintf(fid1,  '%g %g\n',k2d(i),e1plot(h,i));
            fprintf(fid2,  '%g %g\n',k2d(i),e2plot(h,i));
            fprintf(fid3,  '%g %g\n',k2d(i),e3plot(h,i));
            fprintf(fid12, '%g %g\n',k2d(i),e12plot(h,i));
            fprintf(fid123,'%g %g\n',k2d(i),e123plot(h,i));
        end
        fclose('all');
    end
end

 
%% Plot the spectra
COLOR_LINE_ORDER_rainbow_8
if (print_on)
    if (exist(plot_dir) == 0)
        mkdir(plot_dir);
    end
end

line1x = [kappa_max kappa_max];
line1y = [1.0E0 1.0E-10];
line2x = [kappa_max_23 kappa_max_23];
line2y = [1.0E0 1.0E-10];


figure(1)
hp = loglog(k2d,e1plot);
hxl = xlabel('\kappa (1/m)');
hyl = ylabel('E_u(\kappa)');
axis([1.0E-3 1.0E0 1.0E-6 1.0E-0])
PLOTSTYLE_1;

hold on
[X,Y] = ginput(1);
B = log(Y)+(5.0/3.0)*log(X);
Y = exp(-(5.0/3.0)*log(k2d) + B);
hp = loglog(k2d,Y,'Color',[0.75 0.75 0.75]);
PLOTSTYLE_1;
hold off

if (print_on)
    print('-dpng', '-r150',[plot_dir,'/e1x_.png'])
    print('-depsc','-r150',[plot_dir,'/e1x_.eps'])
end

%

figure(2)
hp = loglog(k2d,e2plot);
hxl = xlabel('\kappa (1/m)');
hyl = ylabel('E_v(\kappa)');
axis([1.0E-3 1.0E0 1.0E-6 1.0E-0])
PLOTSTYLE_1;

hold on
[X,Y] = ginput(1);
B = log(Y)+(5.0/3.0)*log(X);
Y = exp(-(5.0/3.0)*log(k2d) + B);
hp = loglog(k2d,Y,'Color',[0.75 0.75 0.75]);
PLOTSTYLE_1;
hold off

if (print_on)
    print('-dpng', '-r150',[plot_dir,'/e2x_.png'])
    print('-depsc','-r150',[plot_dir,'/e2x_.eps'])
end

%

figure(3)
hp = loglog(k2d,e3plot);
hxl = xlabel('\kappa (1/m)');
hyl = ylabel('E_w(\kappa)');
axis([1.0E-3 1.0E0 1.0E-6 1.0E-0])
PLOTSTYLE_1;

hold on
[X,Y] = ginput(1);
B = log(Y)+(5.0/3.0)*log(X);
Y = exp(-(5.0/3.0)*log(k2d) + B);
hp = loglog(k2d,Y,'Color',[0.75 0.75 0.75]);
PLOTSTYLE_1;
hold off

if (print_on)
    print('-dpng', '-r150',[plot_dir,'/e3x_.png'])
    print('-depsc','-r150',[plot_dir,'/e3x_.eps'])
end

%

figure(4)
hp = loglog(k2d,e12plot);
hxl = xlabel('\kappa (1/m)');
hyl = ylabel('E_{uv}(\kappa)');
axis([1.0E-3 1.0E0 1.0E-6 1.0E-0])
PLOTSTYLE_1;

hold on
[X,Y] = ginput(1);
B = log(Y)+(5.0/3.0)*log(X);
Y = exp(-(5.0/3.0)*log(k2d) + B);
hp = loglog(k2d,Y,'Color',[0.75 0.75 0.75]);
PLOTSTYLE_1;
hold off

if (print_on)
    print('-dpng', '-r150',[plot_dir,'/e12x_.png'])
    print('-depsc','-r150',[plot_dir,'/e12x_.eps'])
end

%

figure(5)
hp = loglog(k2d,e123plot);
hxl = xlabel('\kappa (1/m)');
hyl = ylabel('E(\kappa)');
axis([1.0E-3 1.0E0 1.0E-6 1.0E-0])
PLOTSTYLE_1;

hold on
[X,Y] = ginput(1);
B = log(Y)+(5.0/3.0)*log(X);
Y = exp(-(5.0/3.0)*log(k2d) + B);
hp = loglog(k2d,Y,'Color',[0.75 0.75 0.75]);
PLOTSTYLE_1;
hold off

if (print_on)
    print('-dpng', '-r150',[plot_dir,'/e123x_.png'])
    print('-depsc','-r150',[plot_dir,'/e123x_.eps'])
end
