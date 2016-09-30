% BUILD_HSPACE_FROM_CELLS: construct a hierarchical mesh and space, in the
%   parametric domain, from a given set of active cells.
%
%  [hmsh, hspace] = build_hspace_from_cells (dim, degree, num_el0, cells, space_type, truncated)
%
% INPUT:
%
%   dim:        dimension of the parametric domain
%   degree:     degree of the hierarchical splines, it is taken the same in each direction
%   num_el0:    number of elements of level 0, it is taken the same in each direction
%   cells:      cell array with the set of active cells in each level
%   space_type: either 'simplified' (only children functions are activated) or 'standard' (full hierarchical basis)
%   truncated:  a flag to tell whether the basis is truncated or not
%
% OUTPUT:
%
%   hmsh:   object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the coarse space of hierarchical splines (see hierarchical_space)
%
%   The degree and the number of elements are scalars, and are the same in
%    all directions.
%   Regularity is taken to be the maximum (C^{p-1}).
%   The refinement between levels is assumed to be dyadic.
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
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
%

function [hmsh, hspace] = build_hspace_from_cells (dim, degree, initial_num_el, cells, space_type, truncated)

if (nargin < 5 || isempty (space_type))
  space_type = 'standard';
end
if (nargin < 6 || isempty (truncated))
  truncated = false;
end

switch dim
  case 1, problem_data.geo_name = nrbline ([0 0], [1 0]);
  case 2, problem_data.geo_name = 'geo_square.txt';
  case 3, problem_data.geo_name = 'geo_cube.txt';
end


% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = degree*ones(1,dim);         % Degree of the splines
method_data.regularity  = method_data.degree-1;       % Regularity of the splines
method_data.nsub_coarse = initial_num_el*ones(1,dim); % Number of subdivisions
method_data.nsub_refine = 2*ones(1,dim);
method_data.nquad       = method_data.degree+1;       % Points for the Gaussian quadrature rule
method_data.space_type  = space_type;                 % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = truncated;

[hmsh, hspace] = adaptivity_initialize_laplace (problem_data, method_data);

nlevels = numel(cells);

adaptivity_data.flag = 'elements';

for ref = 1:nlevels
    marked = cell (ref,1);
    marked{ref} = cells{ref};
    
    [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);    
end

end
