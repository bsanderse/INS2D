function BC = BC_av_stag3(Nt,Nin,Nb,BC1,BC2,h1,h2)
% total solution u is written as u = Bb*ub + Bin*uin
% the boundary conditions can be written as Bbc*u=ybc
% then u can be written entirely in terms of uin and ybc as:
% u = (Bin-Btemp*Bbc*Bin)*uin + Btemp*ybc, where
% Btemp = Bb*(Bbc*Bb)^(-1)
% Bb, Bin and Bbc depend on type of BC (Neumann/Dirichlet/periodic)


% val1 and val2 can be scalars or vectors with either the value or the
% derivative

% (ghost) points on staggered locations (pressure lines)

% some input checking:
if (Nt~=Nin+Nb)
    error('Number of inner points plus boundary points is not equal to total points');
end

if (Nb==0) % no boundary points, so simply diagonal matrix without boundary contribution
    B1D   = speye(Nt);
    Btemp = spalloc(Nt,2,0);
    ybc1  = zeros(2,1);
    ybc2  = zeros(2,1);
    BC.B1D   = B1D;
    BC.Btemp = Btemp;
    BC.ybc1  = ybc1;
    BC.ybc2  = ybc2;    
    return
end

% boundary conditions
Bbc         = spalloc(Nb,Nt,4);
ybc1_1D     = zeros(Nb,1);
ybc2_1D     = zeros(Nb,1);


if (Nb==6) % normal situation, 2 boundary points

    % boundary matrices
    Bin         = spdiags(ones(Nt,1),-3,Nt,Nin);
    Bb          = spalloc(Nt,Nb,2);
    Bb(1:3,1:3)              = speye(3);
    Bb(end-2:end,end-2:end)  = speye(3);
    
    switch BC1
        case {'dir'}
%             Bbc(1,1)    = 1/2;      % Dirichlet
%             Bbc(1,6)    = 1/2; 
%             Bbc(2,2)    = 1/2;      % Dirichlet uLo
%             Bbc(2,5)    = 1/2;
%             Bbc(3,3)    = 1/2;      % Dirichlet uLo
%             Bbc(3,4)    = 1/2; 
%             ybc1_1D(1)  = 1;        
%             ybc1_1D(2)  = 1;
%             ybc1_1D(3)  = 1;
            Bbc(1,1)    = -1;      
%             Bbc(1,3)    = -1/2;
%             Bbc(1,4)    = 1/2;
            Bbc(1,6)    = 1; 
            Bbc(2,2)    = 1/2;      % Dirichlet uLo
            Bbc(2,5)    = 1/2;
            Bbc(3,3)    = -1;      
            Bbc(3,4)    = 1; 
            ybc1_1D(1)  = 0;        
            ybc1_1D(2)  = 1;
            ybc1_1D(3)  = 0;
        case {'sym'}
%             Bbc(1,1)    = -1;
%             Bbc(1,2)    = 1;
%             ybc1_1D(1)  = h1;     % duLo
            error('not implemented');

        case {'per'}
            Bbc(1:3,1:3)         = -speye(3);
            Bbc(1:3,end-5:end-3) = speye(3);
            Bbc(end-2:end,4:6)   = -speye(3);
            Bbc(end-2:end,end-2:end)  = speye(3);       
        otherwise
            error('not implemented');       
    end

    switch BC2
        case {'dir'}
% standard Dirichlet BC
%             Bbc(end,end)    = 1/2;
%             Bbc(end,end-5)  = 1/2;
%             Bbc(end-1,end-1) = 1/2;
%             Bbc(end-1,end-4) = 1/2;
%             Bbc(end-2,end-2) = 1/2;
%             Bbc(end-2,end-3) = 1/2;
%             ybc2_1D(end-2)   = 1;
%             ybc2_1D(end-1)   = 1;
%             ybc2_1D(end)     = 1;

            Bbc(end,end)     = -1;
            Bbc(end,end-5)   = 1;
            Bbc(end-1,end-1) = 1/2;
            Bbc(end-1,end-4) = 1/2;
            Bbc(end-2,end-2) = -1;
            Bbc(end-2,end-3) = 1;
            ybc2_1D(end-2)   = 0;
            ybc2_1D(end-1)   = 1;
            ybc2_1D(end)     = 0;
%             keyboard
        case {'sym'}
%             Bbc(2,end-1)   = -1;
%             Bbc(2,end)     = 1;
%             ybc2_1D(2)     = h2;     % duUp
            error('not implemented');

        case {'per'} % actually redundant
            Bbc(1:3,1:3)         = -speye(3);
            Bbc(1:3,end-5:end-3) = speye(3);
            Bbc(end-2:end,4:6)   = -speye(3);
            Bbc(end-2:end,end-2:end)  = speye(3);             
        otherwise
            error('not implemented');   
    end


elseif (Nb==1)  % one boundary point 
    
      error('not implemented');

%     Bb          = spalloc(Nt,Nb,2);
%     
%     diagpos = -1;
%     
%     switch BC1
%         case {'dir'}
%             Bbc(1,1)    = 1/2;    
%             Bbc(1,2)    = 1/2;
%             ybc1_1D(1)  = 1;        % uLe    
%             Bb(1,1)     = 1;
%         case {'sym'}
%             Bbc(1,1)    = -1;
%             Bbc(1,2)    = 1;
%             ybc1_1D(1)  = h1;       % duLo
%         case {'per'}
%             Bbc(1,1)    = -1;
%             Bbc(1,end)  = 1;
%             Bb(1,1)     = 1;            
%         otherwise
%             error('not implemented'); 
%     end
%     
%     switch BC2
%         
%         case {'dir'}
%             Bbc(Nb,end-1) = 1/2;
%             Bbc(Nb,end)   = 1/2;   
%             ybc2_1D(1)    = 1;        % uRi   
%             Bb(end,Nb)    = 1;
%         case {'sym'}
%             Bbc(Nb,end-1) = -1;
%             Bbc(Nb,end)   = 1;
%             ybc2_1D(1)    = h2;        % duUp
%         case {'per'}
%             Bbc(1,1)      = -1;
%             Bbc(1,end)    = 1;
%             Bb(1,1)       = 1;
%         otherwise
%             error('not implemented');         
%     end
%     
%     % boundary matrices
%     Bin         = spdiags(ones(Nt,1),diagpos,Nt,Nin);
%  
end

ybc1 = ybc1_1D;
ybc2 = ybc2_1D;
    
Btemp  = Bb*(Bbc*Bb\speye(Nb));     % = inv(Bbc*Bb)
B1D    = Bin - Btemp*Bbc*Bin;

BC.B1D   = B1D;
BC.Btemp = Btemp;
BC.ybc1  = ybc1;
BC.ybc2  = ybc2;