function [contribution] = tangential_contribution(tangentials_,normals_,NV)

tangentials_ok = tangentials_(tangentials_~=-1);
normals_ok = normals_(tangentials_~=-1);
contribution = sparse(tangentials_ok,normals_ok,1,NV,NV);

end