function v = vBC_(x,y,t,options,dare_vectors)
% boundary conditions for v for actuator

% dare_vectors: 1 if you dare to use vectors as inputs, otherwise
if (nargin<5)
    dare_vectors = 0;
end
  
% coordinate left side domain:
x1 = options.grid.x1;

if dare_vectors
    if length(x) ~= length(y)
        error("x and y have different length")
    else
        BCys = y(abs(x-x1)<1e-10);

        f     = 0.5;
        a     = pi/6;
        %     alpha = a*sin(f*t);
        %     v     = sin(alpha)*ones(length(x)*length(y),1);
        alpha = a*sin(BCys-f*t);
        v  = zeros(length(x),1);
        v(abs(x-x1)<1e-10)     = sin(alpha);
    end
elseif ( (length(x)==1 && abs(x-x1)<1e-10)) % || (length(y)==1 && abs(y-0)<1e-10) )
    f     = 0.5;
    a     = pi/6;
%     alpha = a*sin(f*t);
%     v     = sin(alpha)*ones(length(x)*length(y),1);
    alpha = a*sin(y-f*t);
    v     = sin(alpha).*ones(length(x)*length(y),1);
else    
    v = zeros(length(x)*length(y),1);    
end


end