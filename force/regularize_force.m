function f = regularize_force(xmesh,ymesh,xp,yp,panel,param,fk)
% xp,yp are the coordinates of intersections of body with mesh
% panel gives the panel number of the original body panelling with force fk

% for each finite volume intersected by a line segment find the integrated force

nk  = length(xp);
Nx  = length(xmesh);
Ny  = length(ymesh);
imax = Nx;
jmax = Ny;
f    = zeros(Nx,Ny);

% loop over panels
for i=1:nk-1
    k1 = panel(i);
    k2 = panel(i+1);
    
    % find indices of finite volume
    % points inside the FV
    
        xtemp = 0.5*(xp(i)+xp(i+1));
        ytemp = 0.5*(yp(i)+yp(i+1)); 

        % closest u-point=index FV
        [valx indx] = min(abs(xmesh-xtemp));
        [valy indy] = min(abs(ymesh-ytemp));
        
%         xtemp
%         ytemp
%         xmesh(indx)
%         ymesh(indy)
        fact  = 1;
        % check if minimum is unique
        % note: min() returns first value of minimum
        if (indx<imax && abs(xmesh(indx)-xtemp)==abs(xmesh(indx+1)-xtemp))
            indx = indx:indx+1;
            fact = fact+1;
            % panel aligned with vertical mesh face
            warning('panel aligned with vertical face');
%             keyboard
        end
        if (indy<jmax && abs(ymesh(indy)-ytemp)==abs(ymesh(indy+1)-ytemp))
            indy = indy:indy+1;  
            fact = fact+1;
            % panel aligned with horizontal mesh face
            warning('panel aligned with horizontal face');
        end
        
        
    % find integrated force
    if (k2-k1==0) % same panel
      f(indx,indy) = f(indx,indy) + (1/fact)* ( fk(k1)*(param(i+1)-param(i)) );
    elseif (k2-k1==1) % two panels
      f(indx,indy) = f(indx,indy) + (1/fact)* ( fk(k1)*(1-param(i)) + fk(k2)*param(i+1) );
    elseif (k2-k1>1) % more than two panels
      f(indx,indy) = f(indx,indy) + (1/fact)* ( fk(k1)*(1-param(i)) + fk(k2)*param(i+1) + ...
                         sum(fk(k1+1:k2-1)) );
    else
      disp('panels not ordered correctly');    
    end
% f(indx,indy)
% keyboard;
end