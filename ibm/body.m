airfoil_type = '2412';
aoa          = 4;   % in angle of attack
Re_c         = 5e2; % Reynolds number or 0 for inviscid
xfoil        = 1;
rotated      = 1;
% load airfoil properties, Cp and Cf from xfoil or own panel method
% nk panel endpoints, nk-1 panels
load_Cp_Cf;

if (rotated==1)
    [x_k,y_k] = rotate_body(x_k,y_k,aoa);
end

% set coordinates of leading edge (see parameters.m)
% x_c       = L_x/3;
% y_c       = L_y/2;

%% circle:
% R   = 1/2;
% nk  = 51;
% theta = linspace(0,2*pi,nk);
% % theta = theta(1:end-1);
% 
% x_k  = R*cos(theta);
% y_k  = R*sin(theta);

%% square
% nk    = 25;
% L     = 0.5;
% xup   = linspace(L,0,nk)';
% xle   = zeros(nk,1);
% xlo   = linspace(0,L,nk)';
% xri   = L*ones(nk,1);
% x_k   = [xup(1:end-1);xle(1:end-1);xlo(1:end-1);xri(1:end)];
% 
% yup   = L*ones(nk,1);
% yle   = linspace(L,0,nk)';
% ylo   = zeros(nk,1);
% yri   = linspace(0,L,nk)';
% y_k   = [yup(1:end-1);yle(1:end-1);ylo(1:end-1);yri(1:end)];


%% airfoil
%     nk        = 100; 
%     nk2       = floor(nk/2);
%     i         = 1:nk2+1;
%     
%     % clockwise numbering from LE to TE and back
%     x_kb      = sin(0.5*pi*(i-1)*(1/nk2)).^2;   
%     x_ko      = sin(0.5*pi*(1.0-(i-1.0)*(1/nk2))).^2;  
% 
% %     if (~strcmp(airfoil_type,'0012'))
% %         error('not implemented');
% %     end
%     % naca 0012
%     y_kb      = +0.6*(0.29690*x_kb.^0.5 - 0.12600*x_kb - 0.35160*x_kb.^2 + 0.28430*x_kb.^3 -0.10360*x_kb.^4);
%     y_ko      = -0.6*(0.29690*x_ko.^0.5 - 0.12600*x_ko - 0.35160*x_ko.^2 + 0.28430*x_ko.^3 -0.10360*x_ko.^4);
% 
%     % cylinder
%     % y_kb      =  sqrt(0.25-(x_kb-0.5).^2);
%     % y_ko      = -sqrt(0.25-(x_ko-0.5).^2);
% 
%     x_k       = [x_kb(1:end-1) x_ko(1:end)];
%     y_k       = [y_kb(1:end-1) y_ko(1:end)];


%% wave boundary
% number of waves
% n_waves = 1;
% x_k     = linspace(x1,x2,200);
% y_k     = (1/40)*(1-cos(n_waves*2*pi*x_k/L_x)).^2;
% x_k(end+1) = x_k(1);
% y_k(end+1) = y_k(1);

%% shift

x_k       = x_k + x_c;
y_k       = y_k + y_c;
