% HIERARCHICAL_SPACE: constructor of the class for hierarchical spaces.
%
%    function hspace = hierarchical_space (hmsh, space, [space_type, truncated])
%
% INPUT
%    hmsh:       an object of the class hierarchical_mesh (see hierarchical_mesh)
%    space:      the coarsest space, an object of the class sp_bspline (see sp_bpline)
%    space_type: select which kind of hierarchical space to construct. The options are
%                - 'standard',   the usual hierachical splines space (default value)
%                - 'simplified', a simplified basis, were only children of removed functions are activated\
%    truncated:  decide whether the basis will be truncated or not
%
% OUTPUT:
%    hspace: hierarchical_space object, which contains the following fields and methods
% 
%    FIELD_NAME     TYPE                    DESCRIPTION
%    ncomp          (scalar)                number of components of the space
%    type           (string)                'standard' or 'simplified'
%    ndof           (scalar)                total number of active functions 
%    nlevels        (scalar)                the number of levels
%    space_of_level (1 x nlevels)           tensor product space of each level, with 1d evaluations on the mesh of the same level (see sp_bspline)
%    Proj           (hmsh.nlevels-1 x ndim cell-array) 
%                                           the coefficients relating 1D splines of two consecutive levels
%                                           Proj{l,i} is a matrix of size N_{l+1} x N_l where N_l is the number 
%                                           of univariate functions of level l in the direction l, such that
%                                           a function B_{k,l} = \sum_j c^k_j B_{j,l+1}, and c_j = Proj{l,i}(j,k)
%    ndof_per_level (1 x nlevels array)     number of active functions on each level
%    active        (1 x nlevels cell-array) List of active functions on each level
%    coeff_pou     (ndof x 1)               coefficientes to form the partition of the unity in the hierarchical space
%    deactivated   (1 x nlevels cell-array) List of deactivated functions on each level
%    C             (1 x hmsh.nlevels cell-array) Sparse matrices for changing basis. For each level, represent active functions of previous levels
%                                            as linear combinations of splines (active and inactive) of the current level
%    boundary      (2 x ndim array)         a hierarchical space representing the restriction to the boundary
%
%    METHOD NAME
%    hspace_to_vtk:         export the solution to a VTK file, in a structured grid of points
%    hspace_eval:           evaluate the solution in a Cartesian grid of points
%    hspace_eval_hmsh:      evaluate the solution in the quadrature points of a hierarchical mesh
%    hspace_l2_error:       compute the error in L2 norm
%    hspace_h1_error:       compute the error in H1 norm
%    hspace_refine:         refine the hierarchical space
%    hspace_drchlt_l2_proj: compute the boundary degrees of freedom using the L2-projection
%    hspace_check_partition_of_unity: check whether the computed coefficients for the partition of unity are correct 
%    op_gradu_gradv_hier:   assemble the stiffness matrix
%    op_u_v_hier:           assemble the mass matrix
%    op_f_v_hier:           assemble the right-hand side
%
% For an explanation of the 'standard' and the 'simplified' basis for hierarchical splines, read:
%  A. Buffa, E.M. Garau, New refinable spaces and local approximation
%    estimates for hierarchical splines (2015), arxiv:1507.06534v1
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

function hspace = hierarchical_space (hmsh, space, varargin)

default_values = {'standard', false};
default_values(1:numel(varargin)) = varargin;
[space_type, truncated] = default_values{:};

hspace.ncomp = space.ncomp;
hspace.type = space_type;
hspace.truncated = truncated;

hspace.nlevels = 1;
hspace.ndof = space.ndof;
hspace.ndof_per_level = space.ndof;
hspace.space_of_level = space;

hspace.active{1} = (1:space.ndof)';
hspace.deactivated{1} = [];

hspace.coeff_pou = ones (space.ndof, 1);
hspace.Proj = [];
hspace.C{1} = speye (space.ndof);

hspace.dofs = [];

if (~isempty (hmsh.boundary))
  if (hmsh.ndim > 1)
    for iside = 1:2*hmsh.ndim
      boundary = hierarchical_space (hmsh.boundary(iside), space.boundary(iside), space_type, truncated);
      boundary.dofs = space.boundary(iside).dofs;
      hspace.boundary(iside) = boundary;
    end
  elseif (hmsh.ndim == 1)
    hspace.boundary(1).dofs = space.boundary(1).dofs;
    hspace.boundary(2).dofs = space.boundary(2).dofs;
  end
else
  hspace.boundary = [];
end

hspace = class (hspace, 'hierarchical_space');

end