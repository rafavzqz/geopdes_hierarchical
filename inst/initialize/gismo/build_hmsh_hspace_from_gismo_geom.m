% BUILD_HMSH_HSPACE_FROM_GISMO_GEOM: builds the hierarchical mesh and space
% corresponding to a given gismo geometry.
%
% function [hmsh, hspace, method_data] = build_hmsh_hspace_from_gismo_geom (geometry, method_data)
%
% INPUT
%   geometry: a structure built with geo_load from a gsTHBSpline G+smo object. 
%   method_data : a structure with discretization data that contains at least the fields:
%    - nquad:        number of points for Gaussian quadrature rule
%   and possibly (otherwise, those fields will be determined by geometry and returned)
%    - degree:       degree of the spline functions.
%    - regularity:   continuity of the spline functions.
%    - nsub_refine:  number of subelements to be added at each refinement step (2 for dyadic)
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

if ~isfield(method_data, 'degree')
    method_data.degree = geometry.order - 1;
end
if ~isfield(method_data, 'regularity')
    if size(geometry.regularity, 1) == 1
        method_data.regularity = geometry.regularity;
    else
        method_data.regularity = min (geometry.regularity);
    end
end
if isfield(geometry, 'nsub_refine')
    if isfield(method_data, 'nsub_refine') && (~isequal(method_data.nsub_refine, geometry.nsub_refine))
        warning('method_data.nsub_refine is taken as the one describing the G+smo geometry.')
    end
    method_data.nsub_refine = geometry.nsub_refine;
elseif ~isfield(method_data, 'nsub_refine')
    method_data.nsub_refine = 2 * ones (1, length(geometry.order));
    geometry.nsub_refine = method_data.nsub_refine;
else
    geometry.nsub_refine = method_data.nsub_refine;
end
if ~isfield(method_data, 'space_type')
    method_data.space_type = 'standard';
end
if ~isfield(method_data, 'truncated')
    method_data.truncated = 1;
end
if isfield(method_data, 'nsub_coarse') && ~isequal(method_data.nsub_coarse, ones (1, length(geometry.order)))
    warning('method_data.nsub_coarse is not considered when an hmsh is built from a gsTHBSpline object.')
end

% Non hierarchical mesh and space of level 1
knots_geom = geometry.knots{1};
[knots_method, zeta_method] = kntrefine (knots_geom, zeros (1, length(geometry.order)), ...
    method_data.degree, method_data.regularity);
rule = msh_gauss_nodes (method_data.nquad);

% - Build geometrical mesh of level 1
[qn_geom, qw_geom] = msh_set_quad_nodes (knots_geom, rule);
msh_geom = msh_cartesian (knots_geom, qn_geom, qw_geom, geometry);

% - Build mesh and space for the method
[qn, qw] = msh_set_quad_nodes (zeta_method, rule);
msh = msh_cartesian (zeta_method, qn, qw, geometry);
space = sp_bspline (knots_method, method_data.degree, msh);

% Retrieve active elements at each level
[b1, b2, lev_boxes] = geometry.gismo.basis.getBoxes;
cells = boxes_to_elements(b1, b2, lev_boxes, geometry.nsub_refine, msh_geom.nel_dir);

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

