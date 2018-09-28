function BC = BC_general_stag_diff(Nt,Nin,Nb,BC1,BC2,h1,h2)
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

if (Nb==2) % normal situation, 2 boundary points

    % boundary matrices
    Bin         = spdiags(ones(Nt,1),-1,Nt,Nin);
    Bb          = spalloc(Nt,Nb,2);
    Bb(1,1)     = 1;
    Bb(end,Nb)  = 1;
    
    switch BC1
        case {'dir'}
            % second order bc
%             Bbc(1,1)    = 5/16;
%             Bbc(1,2)    = 15/16;
%             Bbc(1,3)    = -5/16;
%             Bbc(1,4)    = 1/16;
            % first order bc
%             Bbc(1,1)    = 3/8;
%             Bbc(1,2)    = 3/4;
%             Bbc(1,3)    = -1/8;    
            % zeroth order (standard mirror conditions)
            Bbc(1,1)    = 1/2;
            Bbc(1,2)    = 1/2;
            ybc1_1D(1)  = 1;        % uLo
        case {'sym'}
            Bbc(1,1)    = -1;
            Bbc(1,2)    = 1;
            ybc1_1D(1)  = h1;     % duLo
        case {'per'}
            Bbc(1,1)    = -1;
            Bbc(1,end-1)= 1;
            Bbc(2,2)    = -1;
            Bbc(2,end)  = 1;         
        otherwise
            error('not implemented');       
    end

    switch BC2
        case {'dir'}
            % second order bc           
%             Bbc(end,end-3)    = 1/16;
%             Bbc(end,end-2)    = -5/16;
%             Bbc(end,end-1)    = 15/16;
%             Bbc(end,end)      = 5/16;        
            % first order bc            
%             Bbc(end,end-2) = -1/8;
%             Bbc(end,end-1) = 3/4;
%             Bbc(end,end)   = 3/8;
            % zeroth order (standard mirror conditions)
            Bbc(end,end-1) = 1/2;
            Bbc(end,end)   = 1/2;            
            ybc2_1D(2)     = 1;     % uUp
        case {'sym'}
            Bbc(2,end-1)   = -1;
            Bbc(2,end)     = 1;
            ybc2_1D(2)     = h2;     % duUp
        case {'per'}
            Bbc(1,1)    = -1;
            Bbc(1,end-1)= 1;
            Bbc(2,2)    = -1;
            Bbc(2,end)  = 1;                
        otherwise
            error('not implemented');   
    end


elseif (Nb==1)  % one boundary point 

    disp('should be unnecessary');
    
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
 
end

ybc1 = ybc1_1D;
ybc2 = ybc2_1D;
    
Btemp  = Bb*(Bbc*Bb\speye(Nb));     % = inv(Bbc*Bb)
B1D    = Bin - Btemp*Bbc*Bin;

BC.B1D   = B1D;
BC.Btemp = Btemp;
BC.ybc1  = ybc1;
BC.ybc2  = ybc2;