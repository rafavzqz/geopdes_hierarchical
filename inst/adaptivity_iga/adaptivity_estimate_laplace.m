% ADAPTIVITY_ESTIMATE_LAPLACE: Computation of a posteriori error indicators for Laplacian problem, using globally smooth (C^1) hierarchical spaces.
%
% We consider the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^n)
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% USAGE:
%
%   est = adaptivity_estimate_laplace (u, hmsh, hspace, problem_data, adaptivity_data)
%
% INPUT:
%
%   u:            degrees of freedom
%   hmsh:         object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace:       object representing the space of hierarchical splines (see hierarchical_space)
%   problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - c_diff:        diffusion coefficient (epsilon in the equation), assumed to be a smooth (C^1) function
%    - grad_c_diff:   gradient of the diffusion coefficient (equal to zero if not present)
%    - f:             function handle of the source term
%   adaptivity_data: a structure with the data for the adaptivity method. In particular, it contains the fields:
%    - flag:          'elements' or 'functions', depending on the refinement strategy.
%    - C0_est:        multiplicative constant for the error indicators
%
%
% OUTPUT:
%
%   est: computed a posteriori error indicators
%           - When adaptivity_data.flag == 'elements': for an element Q, (Buffa and Giannelli, 2016)
%                          est_Q := C0_est * h_Q * (int_Q |f + div(epsilon(x) grad(U))|^2)^(1/2),
%           where h_Q is the local meshsize and U is the Galerkin solution
%           - When adaptivity_data.flag == 'functions': for a B-spline basis function b, (Buffa and Garau, 2016)
%                          est_b := C0_est * h_b * (a_b * int_{supp b} |f + div(epsilon(x) grad(U))|^2 * b)^(1/2),
%           where h_b is the local meshsize, a_b is the coefficient of b for the partition-of-unity, and U is the Galerkin solution
%
% For multipatch domains, with C^0 continuity, jump terms between patches are also considered
%
%        est_Q := C0_est * (h_Q^2 * int_Q |f + div(epsilon(x) grad(U))|^2 + h_Q int_I [epsilon dU/dn]^2 )^(1/2),
% where I is the intersection of the boundary of the element with the interface, 
% and [*] denotes the jump at the interface
%
%        est_b := C0_est * (h_b^2 a_b * int_{supp b} |f + div(epsilon(x) grad(U))|^2 * b + h_b * a_b * int_I [epsilon dU/dn]^2 * b)^(1/2),
% where I is the intersection of the support of the function with the interface
%  
%
%
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2017 Rafael Vazquez
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


function est = adaptivity_estimate_laplace (u, hmsh, hspace, problem_data, adaptivity_data)

if (isfield (adaptivity_data, 'C0_est'))
  C0_est = adaptivity_data.C0_est;
else
  C0_est = 1;
end

[ders, F] = hspace_eval_hmsh (u, hspace, hmsh, {'gradient', 'laplacian'});
dernum = ders{1};
der2num = ders{2};

x = cell (hmsh.rdim, 1);
for idim = 1:hmsh.rdim;
  x{idim} = reshape (F(idim,:), [], hmsh.nel);
end

aux = 0;
valf = problem_data.f (x{:});
val_c_diff = problem_data.c_diff(x{:});
if (isfield (problem_data, 'grad_c_diff'))
  val_grad_c_diff  = feval (problem_data.grad_c_diff, x{:});
  aux = reshape (sum (val_grad_c_diff .* dernum, 1), size(valf));
end
aux = (valf + val_c_diff.*der2num + aux).^2; % size(aux) = [hmsh.nqn, hmsh.nel], interior residual at quadrature nodes

switch lower (adaptivity_data.flag)
  case 'elements',
    w = []; h = [];
    for ilev = 1:hmsh.nlevels
      if (hmsh.msh_lev{ilev}.nel ~= 0)
        w = cat (2, w, hmsh.msh_lev{ilev}.quad_weights .* hmsh.msh_lev{ilev}.jacdet);
        h = cat (1, h, hmsh.msh_lev{ilev}.element_size(:));
      end
    end
    h = h * sqrt (hmsh.ndim);
        
    est = sum (aux.*w);
    est = h.^2 .* est(:);
    
    % Jump terms, only computed for multipatch geometries
    if (isa (hmsh, 'hierarchical_mesh_mp') && hmsh.npatch > 1)
      jump_est = compute_jump_terms (u, hmsh, hspace, problem_data.c_diff, adaptivity_data.flag);
      est = est + h .* jump_est;
    end

    % Neumann boundary conditions
    if (~isempty (problem_data.nmnn_sides))
      nmnn_est = compute_neumann_terms (u, hmsh, hspace, problem_data, adaptivity_data.flag);
      est = est + h .* nmnn_est;
    end
    est = C0_est * sqrt (est);
    

  case 'functions',
    % Compute the mesh size for each level
    ms = zeros (hmsh.nlevels, 1);
    for ilev = 1:hmsh.nlevels
      if (hmsh.msh_lev{ilev}.nel ~= 0)
        ms(ilev) = max (hmsh.msh_lev{ilev}.element_size);
      end
    end
    ms = ms * sqrt (hmsh.ndim);
    
    Nf = cumsum ([0; hspace.ndof_per_level(:)]);
    dof_level = zeros (hspace.ndof, 1);
    for lev = 1:hspace.nlevels
      dof_level(Nf(lev)+1:Nf(lev+1)) = lev;
    end
    coef = ms(dof_level).^2 .* hspace.coeff_pou(:);

    % Residual terms, a_B * h_B^2 * \int_{supp B} |div(epsilon * grad u) + f|^2 B
    est = zeros (hspace.ndof,1);
    ndofs = 0;
    Ne = cumsum([0; hmsh.nel_per_level(:)]);
    for ilev = 1:hmsh.nlevels
      ndofs = ndofs + hspace.ndof_per_level(ilev);
      if (hmsh.nel_per_level(ilev) > 0)
        ind_e = (Ne(ilev)+1):Ne(ilev+1);
        sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
        b_lev = op_f_v (sp_lev, hmsh.msh_lev{ilev}, aux(:,ind_e));
        dofs = 1:ndofs;
        est(dofs) = est(dofs) + hspace.Csub{ilev}.' * b_lev;
      end
    end
    est = coef .* est;

    % Jump terms, only computed for multipatch geometries
    if (isa (hmsh, 'hierarchical_mesh_mp') && hmsh.npatch > 1)
      coef1 = ms(dof_level) .* hspace.coeff_pou(:);
      jump_est = compute_jump_terms (u, hmsh, hspace, problem_data.c_diff, adaptivity_data.flag);
      est = est + coef1 .* jump_est;
    end
    est = C0_est * sqrt (est);
    
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function est = compute_neumann_terms (u, hmsh, hspace, problem_data, flag)
% Compute the terms for boundary sides with Neumann conditions
  
  if (strcmpi (flag, 'elements'))
    est = zeros (hmsh.nel, 1);
  elseif (strcmpi (flag, 'functions'))
    error ('Not implemented for functions, yet')
    est = zeros (hspace.ndof, 1);
  end
  
  if (~isfield (struct (hmsh), 'npatch')) % Single patch case
    for iside = problem_data.nmnn_sides
      if (hmsh.ndim > 1)
% Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
        gside = @(varargin) problem_data.g(varargin{:},iside);
        dofs = hspace.boundary(iside).dofs;
        rhs(dofs) = rhs(dofs) + op_f_v_hier (hspace.boundary(iside), hmsh.boundary(iside), gside);
      else
%         if (iside == 1)
%           x = hmsh.mesh_of_level(1).breaks{1}(1);
%         else
%           x = hmsh.mesh_of_level(1).breaks{1}(end);
%         end
      end
    end
  else % Multipatch case
  end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function est = compute_jump_terms (u, hmsh, hspace, c_diff, flag)
% Compute the jump terms for multipatch geometries, marking either by elements or by functions.

  if (strcmpi (flag, 'elements'))
    est = zeros (hmsh.nel, 1);
  elseif (strcmpi (flag, 'functions'))
    est = zeros (hspace.ndof, 1);
  end

  interfaces = hspace.space_of_level(1).interfaces;

  for iref = 1:numel(interfaces)
% Generate an auxiliary hierarchical mesh, such that the elements on the
%  interface coincide from each side. This is required for integration.
    [hmsh_aux, interface_elements, interface_active_elements] = hmsh_refined_mesh_for_interface (hmsh, interfaces(iref));
    if (hmsh.nel ~= hmsh_aux.nel)
      hspace_aux = hspace_in_finer_mesh (hspace, hmsh, hmsh_aux);
    else
      hspace_aux = hspace;
    end

% Compute the integral of the jump of the normal derivative at the interface.
    if (strcmpi (flag, 'elements'))
      est_edges = integral_term_by_elements (u, hmsh_aux, hspace_aux, interfaces(iref), interface_elements, c_diff);
      est(interface_active_elements(1,:)) = est(interface_active_elements(1,:)) + est_edges;
      est(interface_active_elements(2,:)) = est(interface_active_elements(2,:)) + est_edges;
    elseif (strcmpi (flag, 'functions'))
      est = est + integral_term_by_functions (u, hmsh_aux, hspace_aux, interfaces(iref), interface_elements, c_diff);
    end
  end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function est = integral_term_by_functions (u, hmsh, hspace, interface, interface_elements, coeff)
% Compute the edge integrals for the estimator by functions
%
% OUTPUT
%    est: an array of size hspace.ndof x 1, with the value for each function of
%
%    int_{I \cap supp B} [A dU/dn]^2 B
%
  est = zeros (hspace.ndof, 1);

  patch(1) = interface.patch1;
  patch(2) = interface.patch2;
  side(1) = interface.side1;
  side(2) = interface.side2;

  for lev = 1:hmsh.nlevels
    if (~isempty (interface_elements{lev}{1}))
      ndof_until_lev = sum (hspace.ndof_per_level(1:lev));
      Nelem = cumsum ([0, hmsh.mesh_of_level(lev).nel_per_patch]);
      u_lev = hspace.Csub{lev} * u(1:ndof_until_lev);

      grad_dot_normal = cell (2, 1);
      for ii = [2 1] % This ordering allows me to keep the variables from the first (master) side
        gnum = hspace.space_of_level(lev).gnum{patch(ii)};
        msh_patch_lev = hmsh.mesh_of_level(lev).msh_patch{patch(ii)};
 
% Set of active elements on the patch that are adjacent to the interface
% The numbering in interface_elements is already ordered to make them automatically coincide
        element_list = get_boundary_indices (side(ii), msh_patch_lev.nel_dir, interface_elements{lev}{ii}-Nelem(patch(ii)));

        msh_side_int = msh_boundary_side_from_interior (msh_patch_lev, side(ii));

        msh_side = msh_eval_boundary_side (msh_patch_lev, side(ii), element_list);
        msh_side_aux = msh_evaluate_element_list (msh_side_int, element_list);
 
        sp_bnd = hspace.space_of_level(lev).sp_patch{patch(ii)}.constructor (msh_side_int);
        spp = sp_evaluate_element_list (sp_bnd, msh_side_aux, 'gradient', true);

        grad = sp_eval_msh (u_lev(gnum), spp, msh_side_aux, 'gradient');
        grad_dot_normal{ii} = reshape (sum (grad .* msh_side.normal, 1), msh_side.nqn, msh_side.nel);

        for idim = 1:hmsh.rdim
          x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
        end
        coeffs = coeff (x{:});
        grad_dot_normal{ii} = grad_dot_normal{ii} .* coeffs;
      end

      % Reorder quadrature points, to consider the relative orientation of the patches
      grad_dot_normal{2} = reorder_quad_points (grad_dot_normal{2}, interface, msh_side.nqn_dir);
      grad_dot_normal = grad_dot_normal{1} + grad_dot_normal{2};

% XXXXX I should use a more local numbering, as in the branch localize_Csub
      b_lev = zeros (hspace.space_of_level(lev).ndof, 1);
      b_lev(gnum) = op_f_v (spp, msh_side, grad_dot_normal.^2);
      est(1:ndof_until_lev) = est(1:ndof_until_lev) + hspace.Csub{lev}.' * b_lev;
    end
  end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function est_edges = integral_term_by_elements (u, hmsh, hspace, interface, interface_elements, coeff)
% Compute the edge integrals for the estimator by elements
%
% OUTPUT
%    est: an array of size nedges x 1, with the value for each edge
%
%    int_{e} [A dU/dn]^2
%
% This information must be then transferred to the adjacent elements

  patch(1) = interface.patch1;
  patch(2) = interface.patch2;
  side(1) = interface.side1;
  side(2) = interface.side2;

  est_edges = [];
  for lev = 1:hmsh.nlevels
    if (~isempty (interface_elements{lev}{1}))
      ndof_until_lev = sum (hspace.ndof_per_level(1:lev));
      Nelem = cumsum ([0, hmsh.mesh_of_level(lev).nel_per_patch]);
      u_lev = hspace.Csub{lev} * u(1:ndof_until_lev);

      grad_dot_normal = cell (2, 1);
      for ii = [2 1] % This ordering allows me to keep the variables from the first (master) side
        gnum = hspace.space_of_level(lev).gnum{patch(ii)};
        msh_patch_lev = hmsh.mesh_of_level(lev).msh_patch{patch(ii)};
 
% Set of active elements on the patch that are adjacent to the interface
% The numbering in interface_elements is already ordered to make them automatically coincide
        element_list = get_boundary_indices (side(ii), msh_patch_lev.nel_dir, interface_elements{lev}{ii}-Nelem(patch(ii)));

        msh_side_int = msh_boundary_side_from_interior (msh_patch_lev, side(ii));

        msh_side = msh_eval_boundary_side (msh_patch_lev, side(ii), element_list);
        msh_side_aux = msh_evaluate_element_list (msh_side_int, element_list);
 
        sp_bnd = hspace.space_of_level(lev).sp_patch{patch(ii)}.constructor (msh_side_int);
        spp = sp_evaluate_element_list (sp_bnd, msh_side_aux, 'gradient', true);

        grad = sp_eval_msh (u_lev(gnum), spp, msh_side_aux, 'gradient');
        grad_dot_normal{ii} = reshape (sum (grad .* msh_side.normal, 1), msh_side.nqn, msh_side.nel);

        for idim = 1:hmsh.rdim
          x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
        end
        coeffs = coeff (x{:});
        grad_dot_normal{ii} = grad_dot_normal{ii} .* coeffs;
      end

      % Reorder quadrature points, to consider the relative orientation of the patches
      grad_dot_normal{2} = reorder_quad_points (grad_dot_normal{2}, interface, msh_side.nqn_dir);
      w = msh_side.quad_weights .* msh_side.jacdet;
      dudn_jump = grad_dot_normal{1} + grad_dot_normal{2};
      
      est_edges = cat (2, est_edges, sum (w.*dudn_jump.^2));
    end
  end
  est_edges = est_edges(:);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function field = reorder_quad_points (field, interface, nqn_dir)
% Reorder quadrature points of adjacent patches, to get a corresponding numbering
% The indices for quadrature points are always in second-last position

  indices = num2cell (repmat (':', 1, numel (size (field))));

  ndim = numel (nqn_dir) + 1;
  nqn = prod (nqn_dir);

  if (ndim == 2)
    if (interface.ornt == -1)
      qpoints = nqn:-1:1;
    else
      qpoints = 1:nqn;  
    end
  elseif (ndim == 3)
    if (interface.flag == 1)
      qpoints = reshape (1:nqn, nqn_dir);
    else % I am using nqn_dir from the opposite side (compare with mp_dg_penalty)
      nqn_dir = fliplr (nqn_dir);
      qpoints = reshape (1:nqn, nqn_dir);
      qpoints = qpoints';
    end
    if (interface.ornt1 == -1)
      qpoints = flipud (qpoints);
    end
    if (interface.ornt2 == -1)
      qpoints = fliplr (qpoints);
    end
    qpoints = qpoints(:)';
  end

  indices{end-1} = qpoints;
  field = field(indices{:});

end
