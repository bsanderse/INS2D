function line_int = line_integral(x,y,nx,ny,ds,type)
% evaluate int ln(x+sqrt(x^2+y^2)) dx 


eps = 1e-14;

k = (ny.*y+nx.*x);
q = (ny.*x-nx.*y);

k = (abs(k)>eps).*k;
q = (abs(q)>eps).*q;



if (type==0)
    % area
    line_int = (x./ny).*(0.5*nx.*x-k);
    line_int(abs(ny)<eps) = 0;
    
elseif (type==1)
    % int int 1/r dx dy
    x = (abs(x)>eps).*x;
    line_int = -log( (y+sqrt(x.^2+y.^2)) .^x) + ...
               -log( (ny.*x-nx.*y + sqrt(x.^2+y.^2) ).^k ) + x;
           
elseif (type==2)
    % int int -y/r^2 dx dy
    line_int = (1/2)*ny.*log( (x.^2+y.^2).^q ) - x + ...
               ny.*k.*atan(q./k); 
    % 'repair' atan(0/0)       
    line_int(q==k & abs(q)<eps) = 0 ;
           
elseif (type==3)
    % int int x/r^2 dx dy
%     x*arctan(y/x) - b*( ln(b/x)-(1/2)*ln(1+(y/x)^2) )/(1+a^2) - ...
%           b*a*( arctan((x+a*y)/b)+arctan(1/a) )/(1+a^2)
    x = (abs(x)<eps).*eps + (abs(x)>=eps).*x;
    y = (abs(y)<eps).*eps + (abs(y)>=eps).*y;
    line_int = -x.*atan(y./x) + ...
              + ny.*( log( (k./(ny.*x)).^k ) - (1/2)*(log( (x.^2 +y.^2).^k) - log( (x.^2).^k) ) ) ...
              - nx.*k.*( atan(q./k) - atan(-ny./nx) );   
    % 'repair' atan(0/0)       
    line_int(q==k & abs(q)<eps) = 0;          
    line_int(y==x & abs(y)<eps) = 0;
    % repair ny=0
    line_int(abs(ny)<eps) = 0;    

end
% 

line_int(abs(ds)<eps) = 0;
if (max(isnan(line_int))==1)
    warning('NaN in line integral');
   keyboard;
end
% if (min(isreal(line_int))==0)
%     warning('imaginary part in line integral');
%     keyboard;
% end
% line_int(isnan(line_int)) = 0;
% line_int(isinf(line_int)) = 0;
end