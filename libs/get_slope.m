function slope = get_slope(x,y,n1,n2)
% determine slope (=order) of line in log-log sense
% n1 and n2 (optional) determine interval over which slope is determined; 
% default is n1 = end-1, n2 = end

x = x(:);
y = y(:);

M = length(x);
N = length(y);

stop = 0;

if (M~=N)
    warning('x and y should have same length');
    stop = 1;
end

if (M==1)
    warning('x and y should have at least length 2');
    stop = 1;    
end

if (nargin==2)
    n1 = N-1;
    n2 = N;
end

if (nargin==3)
    if (n1>M)
        error('too many points requested');
    end
    % take n1 to be number of points for which we want slope
    slope = (log(y(2:n1))-log(y(1:n1-1)))./(log(x(2:n1))-log(x(1:n1-1)));
    return
%     warning('wrong number of arguments');
%     stop = 1;    
end

if (stop==0)    
    slope = (log(y(n2))-log(y(n1)))/(log(x(n2))-log(x(n1)));
else
    slope = 0;
end
