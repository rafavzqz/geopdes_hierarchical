% BUILD_GISMO_GEOM_FROM_GEOPDES: builds a hierarchical G+smo geometry
% object gsTHBSpline from a hierarchical GeoPDEs mesh and space, and given
% control points.
%
% function thbs = build_gismo_geom_from_geopdes(hmsh, hspace, pts, [gismo=false])
%
% INPUT:
%   hmsh: hierarchical_mesh object (see hierarchical_mesh.m)
%   hspace: hierarchical_space object (see hierarchical_space.m)
%   pts: array of control points (see argument gismo for the array dimension)
%   gismo: true  if pts are already on the G+smo format of control points (ndof x ndim);
%          false if it is in the GeoPDEs format (ndim+1 x ndof(1) x ... x ndof(ndim))
%          Default: false.
%
% OUTPUT:
%   thbs: gsTHBSpline object (see G+smo documentation)
%   thbsbasis: gsTHBSplineBasis object (see G+smo documentation), the
%       corresponding G+smo hierarchical basis.
%
% Copyright (C) 2019 Ondine Chanon

function [thbs, thbsbasis] = build_gismo_geom_from_geopdes (hmsh, hspace, pts, gismo)

thbsbasis = build_gismo_basis_from_hspace_hmsh(hspace, hmsh);

if nargin < 4 || ~gismo
    if pts(end,:,:,:) ~= ones(size(pts(end,:,:,:)))
        warning('All nurbs weights are put to 1 to build gismo geometry')
    end
    pts = reshape (pts(1:thbsbasis.parDim,:,:,:), thbsbasis.parDim, [])';
end

thbs = gsTHBSpline(thbsbasis, pts, thbsbasis.parDim);

end

