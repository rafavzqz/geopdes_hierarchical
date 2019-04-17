% ADAPTIVITY_INITIALIZE_TH: initialize a hierarchical mesh and a Taylor-Hood pair of spaces with one-single level.
%
% [hmsh, hspace, hspace_press, geometry] = adaptivity_initialize_TH (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - geo_name:     name of the file containing the geometry
%
%  method_data : a structure with discretization data. For this function, it contains the fields:
%    - degree:      degree of the spline functions (pressure space).
%    - regularity:  continuity of the spline functions (pressure space).
%    - nsub_coarse: number of subelements with respect to the geometry mesh (1 leaves the mesh unchanged)
%    - nsub_refine: number of subelements to be added at each refinement step (2 for dyadic)
%    - nquad:       number of points for Gaussian quadrature rule
%    - space_type:  'simplified' (only children of removed functions) or 'standard' (full hierarchical basis)
%    - truncated:   false (classical basis) or true (truncated basis)
%
% OUTPUT:
%    hmsh:         object representing the hierarchical mesh (see hierarchical_mesh)
%    hspace:       object representing a vector-valued space of hierarchical splines (see hierarchical_space)
%    hspace_press: object representing a scalar-valued space of hierarchical splines (see hierarchical_space)
%    geometry:     geometry structure (see geo_load)
%
% The degree of the vector-valued splines is the degree of the scalar-valued splines plus one.
% The regularity of both spaces is the same.
%
% Copyright (C) 2015, 2016 Eduardo M. Garau
% Copyright (C) 2015, 2016, 2017 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [hmsh, hspace, hspace_press, geometry] = adaptivity_initialize_TH (problem_data, method_data)

[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (problem_data.geo_name);
npatch = numel (geometry);

degree_p = method_data.degree;
regularity_p = method_data.regularity;
degree_v = method_data.degree + 1;
regularity_v = method_data.regularity;

msh = cell (1, npatch); 
spp  = cell (1, npatch);
spv  = cell (1, npatch);
for iptc = 1:npatch
  
  [knotsp, zeta] = kntrefine (geometry(iptc).nurbs.knots, method_data.nsub_coarse-1, degree_p, regularity_p);
  rule     = msh_gauss_nodes (method_data.nquad);
  [qn, qw] = msh_set_quad_nodes (zeta, rule);
  msh{iptc}   = msh_cartesian (zeta, qn, qw, geometry(iptc));
  spp{iptc} = sp_bspline (knotsp, degree_p, msh{iptc});
  
  knots_v = kntrefine (geometry(iptc).nurbs.knots, method_data.nsub_coarse-1, degree_v, regularity_v);
  scalar_space = sp_bspline (knots_v, degree_v, msh{iptc});
  for idim = 1:msh{iptc}.ndim
    scalar_spaces{idim} = scalar_space;
  end
  spv{iptc} = sp_vector (scalar_spaces, msh{iptc});
  clear scalar_spaces scalar_space
end

regularity_u = repmat ({method_data.regularity}, 1, msh{iptc}.rdim);
if (npatch == 1)
  msh = msh{1};
  spp = spp{1};
  spv = spv{1};
  hmsh         = hierarchical_mesh (msh, method_data.nsub_refine);
  hspace       = hierarchical_space (hmsh, spv, method_data.space_type, method_data.truncated, regularity_u);
  hspace_press = hierarchical_space (hmsh, spp, method_data.space_type, method_data.truncated, method_data.regularity);
else
  msh = msh_multipatch (msh, boundaries);
  spp = sp_multipatch (spp, msh, interfaces, boundary_interfaces);
  spv = sp_multipatch (spv, msh, interfaces, boundary_interfaces);
  hmsh         = hierarchical_mesh_mp (msh, method_data.nsub_refine);
  hspace       = hierarchical_space_mp (hmsh, spv, method_data.space_type, method_data.truncated, regularity_u);
  hspace_press = hierarchical_space_mp (hmsh, spp, method_data.space_type, method_data.truncated, method_data.regularity);
end

end