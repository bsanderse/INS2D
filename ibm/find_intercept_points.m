function [xbi,ybi,xip,yip,panel_bi] = find_intercept_points(interface_vector,xmesh,ymesh,xk,yk)

eps = 1e-14;


x1 = xmesh(1);
x2 = xmesh(end);
y1 = ymesh(1);
y2 = ymesh(end);

nq = size(interface_vector,1);
nk = length(xk);
[xcol,ycol,Sk,vect,vecn] = panel_properties(xk,yk,1,0);

xbi = zeros(nq,1);
ybi = zeros(nq,1);
panel_bi = zeros(nq,1);
xip = zeros(nq,1);
yip = zeros(nq,1);

% figure
% plot(xk,yk,'kx-');
% hold on

% determine shortest distance to boundary
for q = 1:nq
        
    % loop over interface points   
   i     = interface_vector(q,1);
   j     = interface_vector(q,2);

   xp    = xmesh(i);
   yp    = ymesh(j);


   dx = xk(2:nk)-xk(1:nk-1);
   dy = yk(2:nk)-yk(1:nk-1);
   t  = ((xp-xk(1:nk-1)).*dx + (yp-yk(1:nk-1)).*dy) ./ (Sk.^2);              
   xq = xk(1:nk-1) + t.*dx;
   yq = yk(1:nk-1) + t.*dy;

   % distance vector for all boundary elements
   d  = (xq-xp).^2 + (yq-yp).^2;

%    [val indx] = unique((t>=-eps & t<=1+eps).*d);
   if (min(d)==0) 
       warning('ghost point coinciding with boundary');
       xbi(q) = xp;
       ybi(q) = yp;
       xip(q) = xp;
       yip(q) = yp;

   else
       
       % find closest vertex
       d_vertex = (xp-xk(1:nk-1)).^2 + (yp-yk(1:nk-1)).^2;
       [val_vertex ind_vertex] = min(d_vertex);
       
       % take valid intersections, set others to zero
       [val indx]   = sort( (t>=-eps & t<=1+eps).*d );
       % remove the zeros
       indx2        = find(val);
       % take the minimum by sorting
%        [val2 indx3] = min(d(indx(indx2)));
       [val2 indx3] = sort(d(indx(indx2)));
           
       
       % check if minimum is unique
       if (length(val)>1 && length(indx2)>1)
%            if (abs(val2(indx2(1))-val(indx2(2)))<eps)
           if (abs(val2(1)-val2(2))<eps)               
            disp('minimum not unique');
            % we need an extra ghost cell, which is not 'physical'?               
               k = indx(indx2(1:2));

               % plot the two boundary intercept points
               plot(xq(k),yq(k),'bo-')
               hold on
                              
           end
       end

       
%        plot(xp,yp,'bx')
       % check if distance to boundary is smaller than closest vertex
       if (val2(1) < val_vertex)       

           k = indx(indx2(indx3(1)));
           
           % boundary intercept point
           xbi(q) = xq(k);
           ybi(q) = yq(k);

           % panel associated with boundary intercept
           panel_bi(q) = k;

           % image point
           xip(q) = (xbi(q)-xp)*2 + xp;
           yip(q) = (ybi(q)-yp)*2 + yp;       
%            plot(xbi(q),ybi(q),'ro');
%            plot(xip(q),yip(q),'rx');
       else
           % boundary concave, set boundary intercept equal to closest
           % vertex
           disp('line segment intersection further than closest vertex');           
           % first value in indx3 corresponds to minimum
%            k = indx(indx2(indx3(1)));

           % boundary intercept point
           xbi(q) = xk(ind_vertex);
           ybi(q) = yk(ind_vertex);

           % panel associated with boundary intercept
           panel_bi(q) = ind_vertex; %????????

           % image point
           xip(q) = (xbi(q)-xp)*2 + xp;
           yip(q) = (ybi(q)-yp)*2 + yp;        

%            plot(xbi(q),ybi(q),'go');
%            plot(xip(q),yip(q),'gx');       
       end
       
       % now check if image point not outside domain
       if (xip(q)<x1 || xip(q)>x2 || yip(q)<y1 || yip(q)>y2)
           disp('image point outside domain, trying other point...');
%            keyboard
           if length(indx3)==1
               error('no other point available?');
           end
           k = indx(indx2(indx3(2)));
           % boundary intercept point
           xbi(q) = xq(k);
           ybi(q) = yq(k);

           % panel associated with boundary intercept
           panel_bi(q) = k;

           % image point
           xip(q) = (xq(k)-xp)*2 + xp;
           yip(q) = (yq(k)-yp)*2 + yp;           
       end
   end

%    plot([xp,xq(k)],[yp,yq(k)],'bo-');
%    hold on
%    plot([xp,xip(q)],[yp,yip(q)],'gx-');
               
%    keyboard;
end