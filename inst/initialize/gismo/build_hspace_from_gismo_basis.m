% BUILD_HSPACE_FROM_GISMO_BASIS: builds the GeoPDEs hierarchical space 
% corresponding to a given G+smo hierarchical basis object gsTHBSplineBasis.
%
% function hspace = build_hspace_from_gismo_basis( thbsbasis )
%
% INPUT
%   thbsbasis: gsTHBSplineBasis object (see G+smo documentation)
%
% OUTPUT
%   hspace: hierarchical_space object (see hierarchical_space.m)
%
% Copyright (C) 2019 Ondine Chanon

function hspace = build_hspace_from_gismo_basis( thbsbasis )

coefs = zeros(thbsbasis.size, thbsbasis.parDim);
thbs = gsTHBSpline(thbsbasis, coefs, thbsbasis.parDim);
geometry = geo_load(thbs);
method_data.nquad = 2*(geometry.order-1);

[~, hspace] = build_hmsh_hspace_from_gismo_geom (geometry, method_data);

end

