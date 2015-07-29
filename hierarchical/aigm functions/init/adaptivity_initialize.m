% ADAPTIVITY_INITIALIZE: intialize a hierarchical mesh and a hierarchical
%  space with one-single level.
%
% [hmsh, hspace] = adaptivity_initialize (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - geo_name:     name of the file containing the geometry
%
%  method_data : a structure with discretization data. For this function, it contains the fields:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%    - space_type: 0 (simplified basis), 1 (full basis)
%
% OUTPUT:
%    hmsh:    XXXXXXXXXXXXXXXXXXXXXXXXXXXXxxXXXXXXXXXXXXXXXXXXXX
%    hspace:  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
%    geometry: geometry structure (see geo_load)
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
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

function [hmsh, hspace, geometry] = adaptivity_initialize (problem_data, method_data)

% data_names = fieldnames (problem_data);
% for iopt  = 1:numel (data_names)
%   eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
% end
% data_names = fieldnames (method_data);
% for iopt  = 1:numel (data_names)
%   eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
% end

geometry  = geo_load (problem_data.geo_name);
[knots, zeta] = kntrefine (geometry.nurbs.knots, method_data.nsub-1, method_data.degree, method_data.regularity);
  
rule     = msh_gauss_nodes (method_data.nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);
space    = sp_bspline (knots, method_data.degree, msh);


hmsh     = hierarchical_mesh (msh, geometry);
hspace   = hierarchical_space (hmsh, space, method_data.space_type);

end