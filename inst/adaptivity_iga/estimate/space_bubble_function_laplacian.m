% SPACE_BUBBLE_FUNCTION_LAPLACIAN: computation of the space struct for bubble function estimator
%
% USAGE:
%
%   space_bubble = space_bubble_function_laplacian (hmsh, degree)
%
% INPUT:
%
%   u:            degrees of freedom
%   hmsh:         object representing the hierarchical mesh (see hierarchical_mesh)
%
% OUTPUT:
%
%   space_bubble: computed space struct with bubble functions on active elements
%
%   For more information on the bubble estimators, see Coradello, Antolin and Buffa, CMAME, 2020.
%         
%
% Copyright (C) 2018-2022 Luca Coradello, Rafael Vazquez
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

function sp = space_bubble_function_laplacian (hmsh, degree, varargin)
  
  if numel(varargin) == 0
      iside = -1;
  else
      iside = varargin{1};
  end

  sp.transform = 'grad-preserving';

  sp.space_type = 'spline';
  sp.degree = degree + 1;
  sp.nsh_max = 1;
  for iDim = 1:hmsh.ndim
    sp.nsh_max = sp.nsh_max * (sp.degree(iDim)-1);
  end
  % (TODO: NOT VALID FOR BOUNDARY FUNCTIONS)
  sp.nsh = sp.nsh_max * ones(1,hmsh.nel);
  sp.ndof = hmsh.nel * sp.nsh_max;
  sp.ndof_per_level = [hmsh.nel_per_level] * sp.nsh_max;
%   sp.space_of_level = [];

% Build bubble functions as Bernstein polynomials on a reference element
  ndim = hmsh.ndim;
  for iDim = 1:ndim;
    breaks_reference{iDim} = [0 1];
    knots_reference{iDim} = [zeros(1,sp.degree(iDim)+1) ones(1,sp.degree(iDim)+1)];
  end
  if (ndim == 1)
    geometry = geo_load (nrbline ([0 0], [1 0]));
  elseif (ndim == 2)
    geometry = geo_load (nrb4surf ([0 0], [1 0], [0 1], [1 1]));
  elseif (ndim == 3)
    geometry = geo_load (nrbextrude (nrb4surf ([0 0], [1 0], [0 1], [1 1]), [0 0 1]));
  end
  nqn_dir = hmsh.mesh_of_level(1).nqn_dir;
  rule     = msh_gauss_nodes (nqn_dir);
  [qn, qw] = msh_set_quad_nodes (breaks_reference, rule);
  msh_reference = msh_cartesian (breaks_reference, qn, qw, geometry);
  space_reference = sp_bspline (knots_reference, sp.degree, msh_reference);

  msh_one_elem = msh_evaluate_element_list (msh_reference, 1);
  sp_one_elem = sp_evaluate_element_list_param (space_reference, msh_one_elem, 'value', true, 'gradient', true, 'hessian', false);

  % (TODO: NOT VALID FOR BOUNDARY FUNCTIONS)
  ind_bubbles = cell (ndim, 1);
  for iDim = 1:ndim
    ind_bubbles{iDim} = 2:sp.degree(iDim);
  end
  [ind_bubbles{:}] = ndgrid(ind_bubbles{:});
  bubbles = sub2ind (sp_one_elem.ndof_dir, ind_bubbles{:});
  
  sp_one_elem.shape_functions = sp_one_elem.shape_functions(:,bubbles(:),:);
  sp_one_elem.shape_function_gradients = sp_one_elem.shape_function_gradients(:,:,bubbles(:),:);
%  sp_one_elem.shape_function_hessians = sp_one_elem.shape_function_hessians(:,:,:,bubbles(:),:);
  sp_one_elem.connectivity = sp_one_elem.connectivity(bubbles(:),:);
  
  for iLevel = 1:hmsh.nlevels
    if (hmsh.nel_per_level(iLevel) > 0)
      nel_of_lev = hmsh.nel_per_level(iLevel);
      
      inds = cell (ndim, 1);
      [inds{:}] = ind2sub (hmsh.mesh_of_level(iLevel).nel_dir, hmsh.active{iLevel});

    % Compute element lengths
      length_univ = cellfun (@diff, hmsh.mesh_of_level(iLevel).breaks, 'UniformOutput', false);
      
      lengths = zeros (ndim, 1, 1, nel_of_lev);
      for iDim = 1:ndim
        lengths(iDim,:,:,:) = length_univ{iDim}(inds{iDim});
      end
      ld_first = reshape (lengths, [ndim, 1, 1, 1, nel_of_lev]);
      ld_second = reshape (lengths, [1, ndim, 1, 1, nel_of_lev]);

    % Build space struct of level (TODO: NOT VALID FOR BOUNDARY FUNCTIONS)
      space_of_level.ndof = nel_of_lev * sp.nsh_max;
      space_of_level.ncomp = 1;
      space_of_level.nsh_max = sp.nsh_max;
      space_of_level.nsh = sp.nsh_max * ones(1, nel_of_lev);
      space_of_level.connectivity = reshape (1:space_of_level.ndof, sp.nsh_max, nel_of_lev);
      
    % Multiply basis function derivatives by the correct coefficient, from element length
      shape_funs = repmat (sp_one_elem.shape_functions, [1 1 nel_of_lev]);
      shape_grads = repmat (sp_one_elem.shape_function_gradients, [1 1 1 nel_of_lev]);
%      shape_hess = repmat (sp_one_elem.shape_function_hessians, [1 1 1 1 nel_of_lev]);
      
      space_of_level.shape_functions = shape_funs;
      space_of_level.shape_function_gradients = shape_grads ./ lengths;
%      space_of_level.shape_function_hessians = shape_hess ./ ld_first ./ ld_second;
      
    % Apply grad preserving transform
      space_of_level = sp_grad_preserving_transform (space_of_level, hmsh.msh_lev{iLevel}, 1, 1, 0);
      sp.space_of_level(iLevel) = space_of_level;
    end
  end

end