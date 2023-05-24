function v = vBC_(x,y,t,options,dare_vectors)
% boundary conditions for v for actuator

% dare_vectors: 1 if you dare to use vectors as inputs, otherwise
if (nargin<5)
    dare_vectors = 0;
end

if dare_vectors
    if length(x) ~= length(y)
        error("x and y have different length")
    else
        v = zeros(length(x),1);
    end
else
    v = zeros(length(x)*length(y),1);
end


end