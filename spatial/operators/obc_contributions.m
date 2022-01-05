function C = obc_contributions(tangentials,normals,C,NV,direction,facet_area)

tangentialsA = tangentials(1:end-1);
tangentialsB = tangentials(2:end);

% normal contributions
C = C + direction*.5 *facet_area*sparse(normals,normals,1,NV,NV);
% tangential contributions
C = C + direction*.25*facet_area*tangential_contribution(tangentialsA,normals,NV);
C = C + direction*.25*facet_area*tangential_contribution(tangentialsB,normals,NV);

end