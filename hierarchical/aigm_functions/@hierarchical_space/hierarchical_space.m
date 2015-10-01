% HIERARCHICAL_SPACE: constructor of the class for hierarchical spaces.
%
%    function hspace = hierarchical_space (hmsh, space, [space_type, truncated])
%
% INPUT
%    hmsh:       an object of the class hierarchical_mesh (see hierarchical_mesh)
%    space:      the coarsest space, an object of the class sp_bspline (see sp_bpline)
%    space_type: select which kind of hierarchical space to construct. The options are
%                - 'standard',   the usual hierachical splines space
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
%    globnum_active (ndof x (ndim+1))       global tensor-product numbering of active functions and their corresponding level
%    ndof_per_level (1 x nlevels array)     number of active functions on each level
%    active        (1 x nlevels cell-array) List of active functions on each level
%    coeff_pou     (ndof x 1)               coefficientes to form the partition of the unity in the hierarchical space
%    deactivated   (1 x nlevels cell-array) List of deactivated functions on each level
%    C             (1 x hmsh.nlevels cell-array) Matrices for changing basis (see compute_matrices_for_changing_basis.m)
%    sp_lev        (hmsh.nlevels x 1 cell-array) sp_lev{ilev} is a structure
%    boundary      (2 x ndim array)         a hierarchical space representing the restriction to the boundary
%
%    METHOD NAME
%    hsp_eval       XXXXXXXXX
%    hsp_to_vtk     XXXXXXXXX
%    hsp_l2_error   XXXXXXXXX
%    hsp_h1_error   XXXXXXXXX
%    hspace_refine  XXXXXXXXX
%    operators      XXXXXXXXX
%    plot?
%    split_basis
%    hsp_drchlt_l2_proj?
%
% For an explanation of the 'standard' and the 'simplified' basis for hierarchical splines, read:
%  A. Buffa, E.M. Garau, New refinable spaces and local approximation estimates for hierarchical splines
%  arxiv:1507.06534v1
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

default_values = {0, false}; % Change to 'simplified' and 'standard' for space_type XXXXX
default_values(1:numel(varargin)) = varargin;
[space_type, truncated] = default_values{:};

hspace.ncomp = space.ncomp;
hspace.type = space_type;
hspace.truncated = truncated;

hspace.nlevels = 1;
hspace.ndof = space.ndof;
hspace.ndof_per_level = space.ndof;
hspace.space_of_level = space;

aux = cell (hmsh.ndim, 1);
[aux{:}] = ind2sub ([space.ndof_dir, 1], 1:space.ndof); % The extra 1 makes it work in any dimension
hspace.globnum_active = [ones(space.ndof,1) cell2mat(aux)']; % XXXXX Will this be removed?
hspace.active{1} = (1:space.ndof)';
hspace.deactivated{1} = [];

hspace.coeff_pou = ones (space.ndof, 1);
hspace.Proj = [];
hspace.C{1} = speye (space.ndof);
hspace.sp_lev{1} = sp_evaluate_element_list (hspace.space_of_level(1), hmsh.msh_lev{1}, 'gradient', true,'hessian', true);

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
