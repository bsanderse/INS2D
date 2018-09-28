function area = polygonweightedArea(poly,type)
% poly should be polynomial given by [x,y], size Nx2

% add first point to end
poly = [poly; poly(1,:)];

% delta=[dx,dy]
delta = poly(2:end,:) - poly(1:end-1,:); 

% size of each side:
ds = sqrt(delta(:,1).^2 + delta(:,2).^2);
% normals on each side:
nx = delta(:,2)./ds;
ny = -delta(:,1)./ds;


area = sum(line_integral(poly(2:end,1),poly(2:end,2),nx,ny,ds,type) - ...
           line_integral(poly(1:end-1,1),poly(1:end-1,2),nx,ny,ds,type) ) ;
% we can have imaginary contributions, which normally cancel when summing
% the two line integrals. but when the angle is near 2*pi, the atan
% function might mess up and prevent cancellation. here we manually delete
% any imaginary parts
area = real(area);

if (isinf(area) || isnan(area))
    keyboard;
end
end