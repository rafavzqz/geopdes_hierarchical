% ADAPTIVITY_BUBBLE_ESTIMATOR_KL_SHELL: computation of a posteriori error indicators for Kirchhoff-Love shells, using bubble functions
%
% USAGE:
%
%   est = adaptivity_bubble_estimator_KL_shell (u, hmsh, hspace, problem_data, [method_data, adaptivity_data])
%
% INPUT:
%
%   u:            degrees of freedom
%   hmsh:         object representing the hierarchical mesh (see hierarchical_mesh_mp)
%   hspace:       object representing the space of hierarchical splines (see hierarchical_space_mp_C1)
%   problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - E_coeff:   Young modulus, given as a function handle
%    - nu_coeff:  Poisson ratio, given as a function handle
%    - thickness: thickness of the shell
%    - f:         load, right-hand side of the problem
%    - rotation_sides: if clamped boundary conditions are used (optional)
%   method_data:  structure with data for discretization. If there are clamped condtions (rotation_sides), it must contain
%    - penalty_coeff: penalty coefficient for Nitsche's method.
%   adaptivity_data: structure with data for adaptivity. For this function, it must contain the field:
%    - C0_est:         multiplicative constant for the error indicators
%
%
% OUTPUT:
%
%   est: computed a posteriori error indicators 
%
% WARNING: the current version of the code does not modify the bubble space
%           depending on the boundary conditions
%
%   For more information on the bubble estimators, see Coradello, Antolin and Buffa, CMAME, 2020.
%         
%
% Copyright (C) 2018-2023 Luca Coradello, Rafael Vazquez
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

function estimator = adaptivity_bubble_estimator_KL_shell (u, hmsh, hspace, problem_data, method_data, adaptivity_data)
  
  if (nargin < 6 || ~isfield (adaptivity_data, 'C0_est'))
    C0_est = 1;
  else
    C0_est = adaptivity_data.C0_est;
  end
  E_coeff = problem_data.E_coeff;
  nu_coeff = problem_data.nu_coeff;
  thickness = problem_data.thickness;

  last_dof = cumsum (hspace.ndof_per_level);

  if (isa (hspace.space_of_level(1), 'sp_scalar'))
    error ('Not implemented yet')
    deg = hspace.space_of_level(1).degree;
    space_bubble = space_bubble_function_bilaplacian (hmsh, deg);
    nqn = hmsh.mesh_of_level(1).nqn;
    is_scalar = true;
  elseif (isa (hspace.space_of_level(1), 'sp_multipatch_C1'))
    deg = hspace.space_of_level(1).sp_patch{1}.degree;
    space_bubble = space_bubble_function_bilaplacian_mp (hmsh, deg, true); % true to compute in the parametric domain
    nqn = hmsh.mesh_of_level(1).msh_patch{1}.nqn;
    is_scalar = true;
  elseif (isa (hspace.space_of_level(1), 'sp_vector'))
    deg = hspace.space_of_level(1).scalar_spaces{1}.degree;
    space_bubble = space_bubble_function_bilaplacian (hmsh, deg, true); % true to compute in the parametric domain
    nqn = hmsh.mesh_of_level(1).nqn;
    is_scalar = false;
  end  
  
  estimator = zeros(hmsh.nel, 1);
  
  shifting_vector = cumsum([0 hmsh.nel_per_level]);

  for ilev = 1:hmsh.nlevels
    if (hmsh.nel_per_level(ilev) > 0)
      x = cell (hmsh.rdim,1);
      for idim = 1:hmsh.rdim
        x{idim} = reshape (hmsh.msh_lev{ilev}.geo_map(idim,:,:), nqn, hmsh.nel_per_level(ilev));
      end
%     spu_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', false, 'gradient', false, 'hessian', true);
      if (is_scalar)
        spu_lev = sp_evaluate_element_list_param__ (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true, 'gradient', true, 'hessian', true);
        spu_lev = change_connectivity_localized_Csub (spu_lev, hspace, ilev);
        spu_lev = sp_scalar2vector_param (spu_lev, hmsh.msh_lev{ilev}, 'value', false, 'gradient', true, 'hessian', true);
      else
        spu_lev = sp_evaluate_element_list_param (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true, 'gradient', true, 'hessian', true);
        spu_lev = change_connectivity_localized_Csub (spu_lev, hspace, ilev);
      end
      spv_lev = space_bubble.space_of_level(ilev);
      spv_lev = sp_scalar2vector_param (spv_lev, hmsh.msh_lev{ilev}, 'value', true, 'gradient', true, 'hessian', true);
      
      if (is_scalar) 
        dofs_to_lev = [];
        for icomp = 1:hmsh.rdim
          dofs_to_lev = union (dofs_to_lev, (icomp-1)*hspace.ndof + (1:last_dof(ilev)));
        end
        u_lev = reshape (u(dofs_to_lev), last_dof(ilev), hmsh.rdim);
        solution_of_level = reshape (hspace.Csub{ilev}*u_lev, [], 1);
      else
        solution_of_level = hspace.Csub{ilev}*u(1:last_dof(ilev));
      end

      K_lev = op_KL_shells (spu_lev, spv_lev, hmsh.msh_lev{ilev}, E_coeff(x{:}), nu_coeff(x{:}), thickness);
      b_lev = op_f_v (spv_lev, hmsh.msh_lev{ilev}, problem_data.f(x{:}));           

      K_err_lev = op_KL_shells (spv_lev, spv_lev, hmsh.msh_lev{ilev}, E_coeff(x{:}), nu_coeff(x{:}), thickness);
      if (isfield (problem_data, 'rotation_sides') && ~isempty(problem_data.rotation_sides))
        % Kmat = Kmat + sp_nitsche_KL_rotation (hspace, hmsh, rotation_sides, E_coeff, nu_coeff, thickness, method_data.penalty_coeff);
        [K_lev_nitsche, K_err_lev_nitsche] = nitsche_rotation_estimator (hspace.space_of_level(ilev), ...
          hmsh.mesh_of_level(ilev), hmsh.msh_lev{ilev}, hmsh.boundary, ilev, spv_lev.ndof, hspace.Csub_row_indices{ilev}, ...
          problem_data.rotation_sides, E_coeff, nu_coeff, thickness, method_data.penalty_coeff);
        K_lev = K_lev + K_lev_nitsche;
        K_err_lev = K_err_lev + K_err_lev_nitsche;
      end
      residual_of_level = b_lev - K_lev * solution_of_level;
      error_of_level = K_err_lev \ residual_of_level;

% Compute the estimator from the matrix, avoiding a loop on the elements
      conn = arrayfun (@(x) spv_lev.connectivity(1:spv_lev.nsh(x), x), 1:hmsh.msh_lev{ilev}.nel, 'UniformOutput', false);
      err_elem = cellfun(@(ind) full (error_of_level(ind).' * K_err_lev(ind,ind) * error_of_level(ind)), conn);
      err_elem = sqrt (err_elem);
      
      estimator((shifting_vector(ilev)+1):shifting_vector(ilev+1)) = err_elem;
    end
  end

  estimator = estimator * C0_est;

end

function [K_lev, K_err_lev] = nitsche_rotation_estimator (sp_lev, mesh_of_level, msh_lev, hmsh_bnd, ...
                                    ilev, nbubbles, Csub_row_indices, bnd_sides, E_coeff, nu_coeff, thickness, penalty_coeff)

% Only valid for fixed degree
  rdim = mesh_of_level.rdim;
  ndof_localized_Csub = numel(Csub_row_indices);
  % ndim = mesh_of_level.ndim;
  % degree = sp_lev.sp_patch{1}.degree(1);
  % nbubbles = (degree-2)^ndim*rdim*msh_lev.nel;
  K_lev = sparse (nbubbles, ndof_localized_Csub*rdim);
  K_err_lev = sparse (nbubbles, nbubbles);

  boundaries = mesh_of_level.boundaries;
  Nbnd = cumsum ([0, boundaries.nsides]);
  
  patch_elem_shifting = cumsum ([0 mesh_of_level.nel_per_patch]);
  for iref = bnd_sides
    iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    patch_bnd_shifting = cumsum ([0 hmsh_bnd.mesh_of_level(ilev).nel_per_patch]);

    for ii = 1:numel(iref_patch_list)
      iptc_bnd = iref_patch_list(ii);
      iptc = boundaries(iref).patches(ii);
      iside = boundaries(iref).faces(ii);
      elems_bnd_patch = patch_bnd_shifting(iptc_bnd)+1:patch_bnd_shifting(iptc_bnd+1);
      [~, ~, elements] = intersect (hmsh_bnd.active{ilev}, elems_bnd_patch);

      if (~isempty (elements))
        msh_patch = mesh_of_level.msh_patch{iptc};
% Elements active on the patch and adjacent to the boundary side
        elems_patch = get_volumetric_indices (iside, msh_patch.nel_dir, elements);
        elems_level = patch_elem_shifting(iptc) + elems_patch;
        [~,~,active_elems] = intersect (elems_level, msh_lev.elem_list);

        msh_side = msh_eval_boundary_side (msh_patch, iside, elements);
        msh_side_from_interior = msh_boundary_side_from_interior (msh_patch, iside);

% This is spu, but with the local numbering (level and patch)
        sp_bnd = sp_lev.sp_patch{iptc}.constructor (msh_side_from_interior);
        msh_side_fi = msh_evaluate_element_list (msh_side_from_interior, elements);
        sp_bnd_param = sp_evaluate_element_list_param (sp_bnd, msh_side_fi, 'value', true, 'gradient', true, 'hessian', true);
        sp_bnd_param = sp_scalar2vector_param (sp_bnd_param, msh_side_fi, 'value', true, 'gradient', true, 'hessian', true);

% Bubble functions on the patch, with local (level and patch) numbering/connectivity
        spv_bnd_param = space_bubble_function_on_level_bilaplacian_mp (msh_side_from_interior, msh_side_fi, sp_bnd.degree, true);
        nsh = spv_bnd_param.nsh_max;
        spv_bnd_param.ndof = msh_lev.nel * nsh;
        conn = (active_elems(:)'-1)*nsh + (1:nsh)';
        spv_bnd_param.connectivity = reshape (conn, nsh, msh_side_fi.nel);
        spv_bnd_param = sp_scalar2vector_param (spv_bnd_param, msh_side_fi, 'value', true, 'gradient', true, 'hessian', true);

        for idim = 1:rdim
          x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
        end

% Evaluating parameters on the boundary
        E_bnd  = reshape(E_coeff  (x{:}), msh_side.nqn, msh_side.nel);
        nu_bnd = reshape(nu_coeff (x{:}), msh_side.nqn, msh_side.nel);

        pen_coeff = penalty_coeff * 2^(-ilev);

        K_side = op_nitsche_KL_boundary (sp_bnd_param, spv_bnd_param, msh_side, msh_side_fi, E_bnd, nu_bnd, thickness, pen_coeff);
        K_err_side = op_nitsche_KL_boundary (spv_bnd_param, spv_bnd_param, msh_side, msh_side_fi, E_bnd, nu_bnd, thickness, pen_coeff);

        % [Cpatch, Cpatch_cols_lev] = sp_compute_Cpatch_vector (sp_lev, iptc, rdim);
        [Cpatch, Cpatch_cols_lev] = sp_compute_Cpatch (sp_lev, iptc);
        [~,Csub_rows,Cpatch_cols] = intersect (Csub_row_indices, Cpatch_cols_lev);
        Caux = Cpatch(:,Cpatch_cols); % * hspace.Csub{ilev}(Csub_rows,:);
        Caux = repmat ({Caux}, 1, rdim);
        Caux = blkdiag (Caux{:});

        col_inds = [];
        for icomp = 1:rdim
          col_inds = union (col_inds, (icomp-1)*ndof_localized_Csub + Csub_rows);
        end
        % K_lev(:,Cpatch_cols_lev) = K_lev(:,Cpatch_cols_lev) + K_side * Caux;
        K_lev(:,col_inds) = K_lev(:,col_inds) + K_side * Caux;
        K_err_lev = K_err_lev + K_err_side;
      end

    end
  end
end

% 
function sp = space_bubble_function_on_level_bilaplacian_mp (msh_side_from_interior, msh_side_fi, degree, param)
% % Auxiliary function to create bubble functions on elements adjacent to the
% %  boundary, necessary for Nitsche's method. It creates functions of just
% %  one level.
% % This part of the code is not clean. It is left here in case it is
% % useful, and for reproducibility of results

  if (nargin < 3 || isempty (param))
    param = false;
  end

  nel = msh_side_fi.nel;
  ndim = msh_side_fi.ndim;
  sp.transform = 'grad-preserving';

  sp.space_type = 'spline';
  sp.degree = degree + 1;
  sp.nsh_max = 1;
  for iDim = 1:ndim
    sp.nsh_max = sp.nsh_max * (sp.degree(iDim)-3);
  end
%   (TODO: NOT VALID FOR BOUNDARY FUNCTIONS)
  sp.nsh = ones(1,nel);
  sp.ndof = nel * sp.nsh_max;

% Build bubble functions as Bernstein polynomials on a reference element
  for iDim = 1:ndim
    breaks_reference{iDim} = [0 1];
    knots_reference{iDim} = [zeros(1,sp.degree(iDim)+1) ones(1,sp.degree(iDim)+1)];
  end
  if (ndim == 1)
    geometry = geo_load (nrbline ([0 0], [1 0]));
  elseif (ndim == 2)
    geometry = geo_load (nrb4surf ([0 0], [1 0], [0 1], [1 1]));
  elseif (ndim == 3)
    geometry = geo_load (nrbextrude (nrb4surf ([0 0], [1 0], [0 1], [1 1])));
  end
% XXXX CHANGE quadrature rule to be on the (correct) boundary
%   nqn_dir = msh_side_from_interior.nqn_dir;
%   rule     = msh_gauss_nodes (nqn_dir);
%   [qn, qw] = msh_set_quad_nodes (breaks_reference, rule);
  qn = msh_side_from_interior.qn;
  qw = msh_side_from_interior.qw;
  msh_reference = msh_cartesian (breaks_reference, qn, qw, geometry);
  space_reference = sp_bspline (knots_reference, sp.degree, msh_reference);

  msh_one_elem = msh_evaluate_element_list (msh_reference, 1);
  sp_one_elem = sp_evaluate_element_list_param (space_reference, msh_one_elem, 'value', true, 'gradient', true, 'hessian', true);

  ind_bubbles = cell (ndim, 1);
  for iDim = 1:ndim
    ind_bubbles{iDim} = 3:sp.degree(iDim)-1;
  end
  [ind_bubbles{:}] = ndgrid(ind_bubbles{:});
  bubbles = sub2ind (sp_one_elem.ndof_dir, ind_bubbles{:});

  sp_one_elem.shape_functions = sp_one_elem.shape_functions(:,bubbles(:),:);
  sp_one_elem.shape_function_gradients = sp_one_elem.shape_function_gradients(:,:,bubbles(:),:);
  sp_one_elem.shape_function_hessians = sp_one_elem.shape_function_hessians(:,:,:,bubbles(:),:);
  sp_one_elem.connectivity = sp_one_elem.connectivity(bubbles(:),:);
  sp_one_elem.nsh = numel (bubbles(:));
  sp_one_elem.nsh_max = numel (bubbles(:));

  % Compute element lengths
  lengths = zeros (ndim, 1, 1, nel);

  inds = cell (ndim, 1);
  [inds{:}] = ind2sub (msh_side_fi.nel_dir, msh_side_fi.elem_list);

  length_univ = cellfun (@diff, msh_side_from_interior.breaks, 'UniformOutput', false);

  for iDim = 1:ndim
    lengths(iDim,:,:,1:nel) = length_univ{iDim}(inds{iDim});
  end

  ld_first = reshape (lengths, [ndim, 1, 1, 1, nel]);
  ld_second = reshape (lengths, [1, ndim, 1, 1, nel]);
  
  % Build space struct for active elements (connectivty will have to be changed afterwards)
  sp.ndof = nel * sp.nsh_max;
  sp.ncomp = 1;
  sp.nsh = sp.nsh_max * ones(1, nel);
  sp.connectivity = reshape (1:sp.ndof, sp.nsh_max, nel);

  % Multiply basis function derivatives by the correct coefficient, from element length
  shape_funs = repmat (sp_one_elem.shape_functions, [1 1 nel]);
  shape_grads = repmat (sp_one_elem.shape_function_gradients, [1 1 1 nel]);
  shape_hess = repmat (sp_one_elem.shape_function_hessians, [1 1 1 1 nel]);

  sp.shape_functions = shape_funs;
  sp.shape_function_gradients = shape_grads ./ lengths;
  sp.shape_function_hessians = shape_hess ./ ld_first ./ ld_second;

  % Apply grad preserving transform
  if (~param)
    sp = sp_grad_preserving_transform (sp, msh_side_fi, 1, 1, 1, 1);
  end
  
end