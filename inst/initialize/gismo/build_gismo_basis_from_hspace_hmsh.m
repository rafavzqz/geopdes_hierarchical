% BUILD_GISMO_BASIS_FROM_HSPACE_HMSH: builds a hierarchical G+smo basis
% object gsTHBSplineBasis from a hierarchical GeoPDEs mesh and space.
%
% function thbsbasis = build_gismo_basis_from_hspace_hmsh(hspace, hmsh)
%
% INPUT:
%   hmsh: hierarchical_mesh object (see hierarchical_mesh.m)
%   hspace: hierarchical_space object (see hierarchical_space.m)
%
% OUTPUT:
%   thbsbasis: gsTHBSplineBasis object (see G+smo documentation)
%
% Copyright (C) 2019 Ondine Chanon

function thbsbasis = build_gismo_basis_from_hspace_hmsh(hspace, hmsh)

knots1 = hspace.space_of_level(1).knots;
thbsbasis = gsTHBSplineBasis(knots1, hmsh.ndim);
boxes = elements_to_boxes(hmsh.active, hmsh);
thbsbasis.refineElements(boxes);

end