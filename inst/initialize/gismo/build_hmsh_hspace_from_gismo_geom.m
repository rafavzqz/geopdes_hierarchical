% BUILD_HMSH_HSPACE_FROM_GISMO_GEOM: builds the hierarchical mesh and space
% corresponding to a given gismo geometry.
%
% function [hmsh, hspace, method_data] = build_hmsh_hspace_from_gismo_geom (geometry, method_data)
%
% INPUT
%   geometry: a structure built with geo_load from a gsTHBSpline G+smo
%               object. 
%   method_data : a structure with discretization data that contains the fields:
%    - degree:       degree of the spline functions.
%    - regularity:   continuity of the spline functions.
%    - nsub_refine:  number of subelements to be added at each refinement step (2 for dyadic)
%    - nquad:        number of points for Gaussian quadrature rule
%    - space_type:   'simplified' (only children of removed functions) or 'standard' (full hierarchical basis)
%    - truncated:    false (classical basis) or true (truncated basis)
% 
% OUTPUT
%   hmsh: hierarchical_mesh object (see hierarchical_mesh.m)
%   hspace: hierarchical_space object (see hierarchical_space.m) with minimal
%       regularity (with respect to all the levels) in each parametric direction
%   method_data : a structure with discretization data that contains the fields:
%    - degree:       degree of the spline functions.
%    - regularity:   continuity of the spline functions.
%    - nsub_refine:  number of subelements to be added at each refinement step (2 for dyadic)
%    - nquad:        number of points for Gaussian quadrature rule
%    - space_type:   'simplified' (only children of removed functions) or 'standard' (full hierarchical basis)
%    - truncated:    false (classical basis) or true (truncated basis)
%
% Copyright (C) 2018 Ondine Chanon

function [hmsh, hspace, method_data] = build_hmsh_hspace_from_gismo_geom (geometry, method_data)

method_data.degree = geometry.order - 1;
if size(geometry.regularity, 1) == 1
    method_data.regularity = geometry.regularity;
else
    method_data.regularity = min (geometry.regularity);
end
method_data.space_type = 'standard';
method_data.truncated = 1;

% Non hierarchical mesh and space of level 1
knots1 = geometry.knots{1};
rule = msh_gauss_nodes (method_data.nquad);
[qn, qw] = msh_set_quad_nodes (knots1, rule);
msh = msh_cartesian (knots1, qn, qw, geometry);
space = sp_bspline (knots1, method_data.degree, msh);

% Retrieve active elements at each level
[b1, b2, lev_boxes] = geometry.gismo.basis.getBoxes;
cells = boxes_to_elements(b1, b2, lev_boxes, method_data.nsub_refine, msh.nel_dir);

% Build the hierarchical mesh and space on level 1 only
hmsh = hierarchical_mesh (msh, method_data.nsub_refine);
hspace = hierarchical_space (hmsh, space, method_data.space_type, ...
    method_data.truncated, method_data.regularity);

% Refine hmsh and hspace given the active cells of each level
adaptivity_data.flag = 'elements';
for lev = 1:numel(cells)
    marked_elements = cell(lev,1);
    marked_elements{lev} = setdiff(hmsh.active{lev},cells{lev});
    [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked_elements, ...
        adaptivity_data); 
end

% Update mappings
hmsh = update_mappings(hmsh, geometry);

end

