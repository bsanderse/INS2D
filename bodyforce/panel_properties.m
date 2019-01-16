function [xcol,ycol,Sj,vect,vecn] = panel_properties(xk,yk,contour,airfoil)
% surface, tangential and normal vectors for line segments making a closed contour
% input: xk, yk, vectors with coordinates of panel end points,
%        contour: closed contour (1) or not; if closed we need
%                 xk(1)=xk(end), yk(1)=yk(end)
%        airfoil: standard CCW convention (0) or tangential vector
%                 always pointing in direction of TE (1)
% CAUTION for airfoils with very high angle of attack?



% normally CCW orientation (body always on left side) corresponds with
% correct outward unit normal

if (nargin<3)
    contour = 0;
    airfoil = 0;
end
if (nargin<4)
    airfoil = 0;
end

nk   = length(xk)-1;

% xcol = zeros(nk,1);
% ycol = zeros(nk,1);
% Sj   = zeros(nk,1);
vect = zeros(nk,2);
vecn = zeros(nk,2);

x1   = xk(1:nk);
x2   = xk(2:nk+1);
y1   = yk(1:nk);
y2   = yk(2:nk+1);

xcol = 0.5*(x1+x2);
ycol = 0.5*(y1+y2);
Sj   = sqrt((x2-x1).^2 + (y2-y1).^2);


if (airfoil==0)
    % normal convention, CCW contour
    vect(:,1)   = (x2-x1)./Sj;   
    vect(:,2)   = (y2-y1)./Sj;
    vecn(:,1)   = vect(:,2);
    vecn(:,2)   = -vect(:,1);    
elseif (airfoil==1)
    
    % have the tangential vector always pointing from LE to TE on both
    % upper and lower side
    
    if (contour == 1)
        % define leading edge as point with minimum x_k; this point belongs
        % to upper side 
        [val ile] = min(xk);
        % there is a change that point before or after has same value:
        if (xk(ile)==xk(ile-1))
            ile = ile-1;
        end
        vect(:,1)   = (x2-x1)./Sj;
        vect(:,2)   = (y1-y2)./Sj;
        vecn(:,1)   = -vect(:,2);
        vecn(:,2)   = -vect(:,1);
        % force pointing towards TE
        vect(:,1)   = abs(vect(:,1)); % x-component always positive 
        vect(ile:end,2) = -vect(ile:end,2); 
    elseif (contour == 0)
        vect(:,1)   = (x2-x1)./Sj;
        vect(:,2)   = (y1-y2)./Sj;
        vecn(:,1)   = -vect(:,2);
        vecn(:,2)   = -vect(:,1);
        vect(:,1)   = abs(vect(:,1));  
    end
end