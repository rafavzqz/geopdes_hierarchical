% HIERARCHICAL_SPACE_MP_C1: constructor of the class for hierarchical C1 multipatch spaces.
%
%    function hspace = hierarchical_space_mp_C1 (hmsh, space, [space_type, truncated, regularity])
%
% INPUT
%    hmsh:       an object of the class hierarchical_mesh_mp (see hierarchical_mesh_mp)
%    space:      the coarsest space, an object of the class sp_multipatch (see sp_multipatch_C1)
%    space_type: select which kind of hierarchical space to construct. The options are
%                - 'standard',   the usual hierachical splines space (default value)
%                - 'simplified', a simplified basis, where only children of removed functions are activated
%    truncated:  decide whether the basis will be truncated or not (not truncated by default)
%    regularity: will be used for refinement. For vectors, it should be given in a cell array. By default it is degree minus two
%
% OUTPUT:
%    hspace: hierarchical_space_mp_C1 object, which contains the following fields and methods
% 
%    FIELD_NAME     TYPE                    DESCRIPTION
%    ncomp          (scalar)                number of components of the space
%    type           (string)                'standard' or 'simplified'
%    truncated      (logical)               truncated or non-truncated basis
%    ndof           (scalar)                total number of active functions 
%    nlevels        (scalar)                the number of levels
%    space_of_level (1 x nlevels)           C1 multipatch space of each level (see sp_multipatch_C1)
%    Proj           (hmsh.nlevels-1 x npatch cell-array) 
%                                           the coefficients relating 1D splines of two consecutive levels for each patch
%                                           Proj{l,i} is a cell-array of dimension ndim, with the information for
%                                           the univariate Projectors on the i-th patch (see also hierarchical_space)
%    Proj0                                  Like Proj, for univariate splines of degree p, regularity r+1
%    Proj1                                  Like Proj, for univariate splines of degree p-1, regularity r
%    ndof_per_level (1 x nlevels array)     number of active functions on each level
%    active        (1 x nlevels cell-array) List of active functions on each level
%    deactivated   (1 x nlevels cell-array) List of deactivated functions on each level
%    Csub          (1 x hmsh.nlevels cell-array) Sparse matrices for changing basis. For each level, represent active functions of previous levels
%                                            as linear combinations of splines (active and inactive) of the current level
%    Csub_row_indices (1 x hmsh.nlevels cell-array) indices of the rows stored in Csub. 
%                                            This allows to save memory space.
%
%    METHOD NAME
%    Methods for post-processing, which require a computed vector of degrees of freedom
%      sp_to_vtk:             export the solution to a VTK file, in a structured grid of points on each patch
%      sp_eval:               evaluate the solution in a Cartesian grid of points on each patch
%      hspace_eval_hmsh:      evaluate the solution in the quadrature points of the corresponding hierarchical mesh
%      sp_l2_error:           compute the error in L2 norm
%      sp_h1_error:           compute the error in H1 norm
%
%    Methods for matrix and vector assembly, and boundary conditions
%      op_gradu_gradv_hier:   assemble the stiffness matrix
%      op_u_v_hier:           assemble the mass matrix
%      op_f_v_hier:           assemble the right-hand side
%      sp_drchlt_l2_proj:     compute the boundary degrees of freedom using the L2-projection
%
%    Methods to obtain a relation between different levels
%      hspace_get_children:   get the children of a given set of functions
%      hspace_get_parents:    get the parents of a given set of functions
%      hspace_subdivision_matrix: compute the matrix of basis change (the same stored in Csub)
%
%    Other methods
%      hspace_refine:         refine the hierarchical space
%      hspace_check_partition_of_unity: check whether the computed coefficients
%                             for the partition of unity are correct (used for debugging)
%      hspace_add_new_level:  add a new level, initialized without active functions
%      hspace_remove_empty_level: remove the finest level, if it is empty
%
% For details about the hierarchical multipatch C1 space:
%    C. Bracco, C. Giannelli, M. Kapl, R. Vazquez, Adaptive isogeometric
%    methods with C1 (truncated) hierarchical splines on planar multi-patch
%    domains, Math. Models Meth. Appl. Sci. (2023)
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2017 Rafael Vazquez
% Copyright (C) 2019-2024 Cesare Bracco, Rafael Vazquez
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

function hspace = hierarchical_space_mp_C1 (hmsh, space, varargin)

if (~isa (space, 'sp_multipatch_C1'))
  error ('The zero level space should be of class sp_multipatch_C1')
end

if (isa (space.sp_patch{1}, 'sp_scalar'))
  is_scalar = true;
  regularity = space.sp_patch{1}.degree - 2;
elseif (isa (space.sp_patch{1}, 'sp_vector'))
  error ('Only done for scalars, so far')
else
  error ('Unknown space type')
end

default_values = {'standard', false, regularity};
default_values(1:numel(varargin)) = varargin;
[space_type, truncated, regularity] = default_values{:};

% This is done to simplify things for the boundary and gluing patches
if (check_regularity (regularity, space, is_scalar))
  error ('For multipatch geometries, we force the regularity to be the same in all directions')
end

hspace.ncomp = space.ncomp;
hspace.type = space_type;
hspace.truncated = truncated;

hspace.nlevels = 1;
hspace.ndof = space.ndof;
hspace.ndof_per_level = space.ndof;
hspace.space_of_level = space;

hspace.active{1} = (1:space.ndof)';
hspace.deactivated{1} = [];

hspace.Proj = cell (0, hmsh.npatch);
hspace.Csub{1} = speye (space.ndof);
hspace.Csub_row_indices{1} = 1:space.ndof;

hspace.Proj0 = cell (0, hmsh.npatch);
hspace.Proj1 = cell (0, hmsh.npatch);

hspace.dofs = [];

hspace.regularity = regularity;
hspace = class (hspace, 'hierarchical_space_mp_C1');

end


function error_flag = check_regularity (regularity, space, is_scalar)
% Check if the regularity is the same in all directions.
% For curl-preserving and div-preserving spaces, they should follow the
%  standard relation given by the De Rham diagram, that is, the first space
%  in the sequence should have the same regularity in all directions

error_flag = false;
if (is_scalar)
  if (~all (regularity == regularity(1)))
    error_flag = true;
  end
else
  if (strcmpi (space.sp_patch{1}.transform, 'grad-preserving'))
    for icomp = 1:space.sp_patch{1}.ncomp
      if (~all (regularity{icomp} == regularity{icomp}(1)))
        error_flag = true;
      end
    end
  elseif (strcmpi (space.sp_patch{1}.transform, 'curl-preserving'))
    for icomp = 1:space.sp_patch{1}.ncomp_param
      reg = regularity{icomp};
      reg(icomp) = reg(icomp) + 1;
      if (~all (reg == reg(1)))
        error_flag = true;
      end
    end
  elseif (strcmpi (space.sp_patch{1}.transform, 'div-preserving'))
    for icomp = 1:space.sp_patch{1}.ncomp_param
      reg = regularity{icomp};
      reg(icomp) = reg(icomp) - 1;
      if (~all (reg == reg(1)))
        error_flag = true;
      end
    end
  end
end

end