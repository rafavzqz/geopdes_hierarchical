% ADAPTIVITY_ESTIMATE_LINEAR_EL_MIXED: Computation of a posteriori error indicators for linear elasticity problem, using globally smooth (C^1) hierarchical spaces,
% in the case of mixed formulation.
%
% We consider the linear elasticity problem
%
%      - div (sigma(u)) = f    in Omega = F((0,1)^n)
%      sigma(u) \cdot n = g    on Gamma_N
%                     u = h    on Gamma_D
%
% solve with a mixed formulation
%
% USAGE:
%
%   est = adaptivity_estimate_linear_el_mixed (u, press, hmsh, hspace, hspace_p, problem_data, adaptivity_data)
%
% INPUT:
%
%   u:            degrees of freedom for the displacement
%   press:        degrees of freedom for the multiplier
%   hmsh:         object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace:       object representing the space of hierarchical splines (see hierarchical_space)
%   hspace_p:     object representing the space of hierarchical splines for the multiplier (see hierarchical_space)
%   problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - lambda_lame, mu_lame: function handles of the Lame parameters
%    - f:             function handle of the source term
%   adaptivity_data: a structure with the data for the adaptivity method. In particular, it contains the fields:
%    - flag:          'elements' or 'functions', depending on the refinement strategy.
%    - C0_est:        multiplicative constant for the error indicators
%
%
% OUTPUT:
%
%   est: computed a posteriori error indicators
%           - (Buffa and Giannelli, 2016) When adaptivity_data.flag == 'elements': for an element Q,
%                          est_Q := C0_est*h_Q*(int_Q |f + div(epsilon(x) grad(U))|^2)^(1/2),
%           where h_Q is the local meshsize and U is the Galerkin solution
%           - (Buffa and Garau, 2016) When adaptivity_data.flag == 'functions': for a B-spline basis function b,
%                          est_b := C0_est*h_b*(int_{supp b} a_b*|f + div(epsilon(x) grad(U))|^2*b)^(1/2),
%           where h_b is the local meshsize, a_b is the coefficient of b for the partition-of-unity, and U is the Galerkin solution
%
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2017, 2018 Cesare Bracco, Rafael Vazquez
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


function est = adaptivity_estimate_linear_el_mixed (u, press, hmsh, hspace_u, hspace_p, problem_data, adaptivity_data)

if (isfield(adaptivity_data, 'C0_est'))
    C0_est = adaptivity_data.C0_est;
else
    C0_est = 1;
end

[est1, est2] = compute_residual_terms (u, press, hmsh, hspace_u, hspace_p, problem_data, adaptivity_data.flag);

if (isa (hmsh, 'hierarchical_mesh_mp') && hmsh.npatch > 1)
  jump_est = compute_jump_terms (u, press, hmsh, hspace_u, hspace_p, problem_data.mu_lame, adaptivity_data.flag);
else
  jump_est = 0;
end

if (~isempty (problem_data.nmnn_sides))
  nmnn_est = compute_neumann_terms (u, press, hmsh, hspace_u, hspace_p, problem_data, adaptivity_data.flag);
else
  nmnn_est = 0;
end

switch lower (adaptivity_data.flag)
  case 'elements'
    h = [];
    for ilev = 1:hmsh.nlevels
      if (hmsh.msh_lev{ilev}.nel ~= 0)
        h = cat (1, h, hmsh.msh_lev{ilev}.element_size(:));
      end
    end
    h = h * sqrt (hmsh.ndim);
    
    est1 = h.^2 .* est1(:) + h .* (jump_est + nmnn_est);
    
  case 'functions'
    % Compute the mesh size for each level
    ms = zeros (hmsh.nlevels, 1);
    for ilev = 1:hmsh.nlevels
      if (hmsh.msh_lev{ilev}.nel ~= 0)
        ms(ilev) = max (hmsh.msh_lev{ilev}.element_size);
      else
        ms(ilev) = 0;
      end
    end
    ms = ms * sqrt (hmsh.ndim);
    
    Nf = cumsum ([0; hspace_u.ndof_per_level(:)]);
    dof_level = zeros (hspace_u.ndof, 1);
    for lev = 1:hspace_u.nlevels
      dof_level(Nf(lev)+1:Nf(lev+1)) = lev;
    end
    coef = ms(dof_level).^2 .* hspace_u.coeff_pou(:);
    coef1 = ms(dof_level) .* hspace_u.coeff_pou(:);
    
    est1 = coef .* est1 + coef1 .* (jump_est + nmnn_est);
    est2 = hspace_p.coeff_pou(:) .* est2;
end
est = C0_est .* sqrt (est1(:) + est2(:));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [est1, est2] = compute_residual_terms (u, press, hmsh, hspace_u, hspace_p, problem_data, flag)

  [ders2_u, F] = hspace_eval_hmsh (u, hspace_u, hmsh, {'hessian', 'divergence'});
  [ders2_u, divergence_u] = deal (ders2_u{:});

  [gradient_p, F] = hspace_eval_hmsh (press, hspace_p, hmsh, {'value', 'gradient'});
  [pressure, gradient_p] = deal (gradient_p{:});

  x = cell (hmsh.rdim, 1);
  for idim = 1:hmsh.rdim;  %rdim is the dimension of the physical domain
    x{idim} = reshape (F(idim,:), [], hmsh.nel);
  end

% Check whether lambda and mu are constants
  if (numel (unique (problem_data.lambda_lame (x{:}))) > 1) || ...
    (numel (unique (problem_data.mu_lame (x{:}))) > 1)
    warning ('We assume that the Lame coefficients are constants')
  end

  lambda = problem_data.lambda_lame(1); %we assume that lambda and mu are constants
  mu = problem_data.mu_lame(1);

  valf = problem_data.f (x{1},x{2});
  for ii = 1:hspace_u.ncomp
    partials_a = 0;
    partials_b = 0;
    for jj = 1:hspace_u.ncomp
      if (jj ~= ii)
        partials_a = partials_a + reshape (ders2_u(jj,ii,jj,:,:), [], hmsh.nel); %mixed derivatives of all components (except h-th component)
        partials_b = partials_b + reshape (ders2_u(ii,jj,jj,:,:), [], hmsh.nel); %second derivatives of h-th component (except w.r.t. h-th variable)
      end
    end
    divergence(ii,:,:) = (2*mu) * reshape(ders2_u(ii,ii,ii,:,:), [], hmsh.nel) + mu*partials_a + mu*partials_b + reshape(gradient_p(ii,:,:), [], hmsh.nel);
  end
  aux1 = (valf + divergence).^2;  %residual 1
  aux2 = ((1/lambda)*pressure - divergence_u).^2;   %residual 2

  Ne = cumsum([0; hmsh.nel_per_level(:)]);
  if (strcmpi (flag, 'elements'))
    for ilev = 1:hmsh.nlevels
      if (hmsh.msh_lev{ilev}.nel ~= 0)
        ind_e = Ne(ilev)+1:Ne(ilev+1);
        w = hmsh.msh_lev{ilev}.quad_weights .* hmsh.msh_lev{ilev}.jacdet;
        aux1_elems = reshape (sum (aux1(:,:,ind_e)), [], numel(ind_e));
        est1(ind_e) = sum (w .* aux1_elems);
        est2(ind_e) = sum (w .* aux2(:,ind_e));
      end
    end

  elseif (strcmpi (flag, 'functions'))
    error ('The estimator for the mixed formulation is not implemented for functions: decide the output')
    
    ndofs = 0;
    for ilev = 1:hmsh.nlevels
      ndofs = ndofs + hspace_u.ndof_per_level(ilev);
      ndofs_p = ndofs_p + hspace_p.ndof_per_level(ilev);
      if (hmsh.nel_per_level(ilev) > 0)
        ind_e = (Ne(ilev)+1):Ne(ilev+1);
        sp_lev = sp_evaluate_element_list (hspace_u.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
        b_lev = op_f_v (sp_lev, hmsh.msh_lev{ilev}, aux1(:,:,ind_e));
        dofs = 1:ndofs;
        est1(dofs) = est1(dofs) + hspace_u.Csub{ilev}.' * b_lev;

        sp_lev = sp_evaluate_element_list (hspace_p.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
        b_lev = op_f_v (sp_lev, hmsh.msh_lev{ilev}, aux2(:,ind_e));
        dofs = 1:ndofs_p;
        est2(dofs) = est2(dofs) + hspace_p.Csub{ilev}.' * b_lev;
      end
    end
  end
end
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function est = compute_neumann_terms (u, press, hmsh, hspace, hspace_p, problem_data, flag)
% Compute the terms for boundary sides with Neumann conditions
  
  if (strcmpi (flag, 'elements'))
    est = zeros (hmsh.nel, 1);
  elseif (strcmpi (flag, 'functions'))
    est = zeros (hspace.ndof, 1);
  end

  if (~isfield (struct (hmsh), 'npatch')) % Single patch case
    for iside = problem_data.nmnn_sides
      gside = @(varargin) problem_data.g(varargin{:},iside);
      hmsh_sfi = hmsh_boundary_side_from_interior (hmsh, iside);
      shifting_indices = cumsum ([0 hmsh.boundary(iside).nel_per_level]);
      vol_shifting_indices = cumsum ([0 hmsh.nel_per_level]);
      last_dof = cumsum (hspace.ndof_per_level);
      last_dofp = cumsum (hspace_p.ndof_per_level);

      ndofs = cumsum (hspace.boundary(iside).ndof_per_level);
      ndofs_p = cumsum (hspace_p.boundary(iside).ndof_per_level);
      for ilev = 1:hmsh.boundary(iside).nlevels
        if (hmsh.boundary(iside).nel_per_level(ilev) > 0)
          elements = shifting_indices(ilev)+1:shifting_indices(ilev+1);
          msh_side = hmsh_eval_boundary_side (hmsh, iside, elements);
          msh_side_from_interior = hmsh_sfi.mesh_of_level(ilev);

          sp_bnd = hspace.space_of_level(ilev).constructor (msh_side_from_interior);
          msh_side_from_interior_struct = msh_evaluate_element_list (msh_side_from_interior, hmsh_sfi.active{ilev});
          sp_bnd_struct = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', false, 'gradient', true);
          sp_bnd = hspace_p.space_of_level(ilev).constructor (msh_side_from_interior);
          spp_bnd_struct = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', true);

          grad = sp_eval_msh (hspace.Csub{ilev}*u(1:last_dof(ilev)), sp_bnd_struct, msh_side_from_interior_struct, 'gradient');
          gradt = permute (grad, [2 1 3 4]);
          normal = reshape (msh_side.normal, 1, [], msh_side.nqn, msh_side.nel);
          eps_normal = reshape (sum (bsxfun (@times, grad + gradt, normal), 2), [], msh_side.nqn, msh_side.nel);
          
          press_val = reshape (sp_eval_msh (hspace_p.Csub{ilev}*press(1:last_dofp(ilev)), spp_bnd_struct, ...
                                    msh_side_from_interior_struct, 'value'), 1, msh_side.nqn, msh_side.nel);
          p_normal = reshape (bsxfun (@times, press_val, msh_side.normal), [], msh_side.nqn, msh_side.nel);

          for idim = 1:hmsh.rdim
            x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
          end
          coeff = reshape (problem_data.mu_lame (x{:}), 1, msh_side.nqn, msh_side.nel);
          coeff = (bsxfun (@times, eps_normal, coeff) + p_normal - gside(x{:})).^2;

          if (strcmpi (flag, 'elements'))
            est_level = sum (reshape (sum (coeff, 1), msh_side.nqn, msh_side.nel));
            inds_level = get_volumetric_indices (iside, hmsh.mesh_of_level(ilev).nel_dir, hmsh_sfi.active{ilev});
            [~,~,inds] = intersect (inds_level, hmsh.active{ilev});
            indices = vol_shifting_indices(ilev) + inds;
            est(indices) = est_level;
%          elseif (strcmpi (flag, 'functions'))
%            msh_side = msh_evaluate_element_list (hmsh.boundary(iside).mesh_of_level(ilev), hmsh.boundary(iside).active{ilev});
%            sp_bnd = sp_evaluate_element_list (hspace.boundary(iside).space_of_level(ilev), msh_side, 'value', true);
%            est_level = op_f_v (sp_bnd, msh_side, coeff);
%            dofs = hspace.boundary(iside).dofs(1:ndofs(ilev));
%            est(dofs) = est(dofs) + hspace.boundary(iside).Csub{ilev}.' * est_level;
          end
        end
      end
    end
  else % Multipatch case
    boundaries = hmsh.mesh_of_level(1).boundaries;
    Nbnd = cumsum ([0, boundaries.nsides]);
    last_dof = cumsum (hspace.ndof_per_level);
    last_dofp = cumsum (hspace_p.ndof_per_level);
    
    ndofs = cumsum (hspace.boundary.ndof_per_level);
    for iref = problem_data.nmnn_sides
      iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
      gside = @(varargin) problem_data.g(varargin{:},iref);
      vol_shifting_indices = cumsum ([0 hmsh.nel_per_level]);
      for ilev = 1:hmsh.boundary.nlevels
        u_lev = hspace.Csub{ilev}*u(1:last_dof(ilev));
        press_lev = hspace_p.Csub{ilev}*press(1:last_dofp(ilev));
        patch_bnd_shifting = cumsum ([0 hmsh.boundary.mesh_of_level(ilev).nel_per_patch]);
        patch_shifting = cumsum ([0 hmsh.mesh_of_level(ilev).nel_per_patch]);

        for ii = 1:numel(iref_patch_list)
          iptc_bnd = iref_patch_list(ii);
          iptc = boundaries(iref).patches(ii);
          iside = boundaries(iref).faces(ii);
          elems_patch = patch_bnd_shifting(iptc_bnd)+1:patch_bnd_shifting(iptc_bnd+1);
          [~, ~, elements] = intersect (hmsh.boundary.active{ilev}, elems_patch);
          
          if (~isempty (elements))
            msh_patch = hmsh.mesh_of_level(ilev).msh_patch{iptc};

            msh_side = msh_eval_boundary_side (msh_patch, iside, elements);
            msh_side_from_interior = msh_boundary_side_from_interior (msh_patch, iside);

            sp_bnd = hspace.space_of_level(ilev).sp_patch{iptc}.constructor (msh_side_from_interior);
            msh_side_from_interior_struct = msh_evaluate_element_list (msh_side_from_interior, elements);
            sp_bnd_struct = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', false, 'gradient', true);
            sp_bnd = hspace_p.space_of_level(ilev).sp_patch{iptc}.constructor (msh_side_from_interior);
            spp_bnd_struct = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', true);

            u_patch = u_lev(hspace.space_of_level(ilev).gnum{iptc});
            press_patch = press_lev(hspace_p.space_of_level(ilev).gnum{iptc});
            
            grad = sp_eval_msh (u_patch, sp_bnd_struct, msh_side_from_interior_struct, 'gradient');
            gradt = permute (grad, [2 1 3 4]);
            normal = reshape (msh_side.normal, 1, [], msh_side.nqn, msh_side.nel);
            eps_normal = reshape (sum (bsxfun (@times, grad + gradt, normal), 2), [], msh_side.nqn, msh_side.nel);
          
            press_val = reshape (sp_eval_msh (press_patch, spp_bnd_struct, ...
                                    msh_side_from_interior_struct, 'value'), 1, msh_side.nqn, msh_side.nel);
            p_normal = reshape (bsxfun (@times, press_val, msh_side.normal), [], msh_side.nqn, msh_side.nel);

            for idim = 1:hmsh.rdim
              x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
            end
            coeff = reshape (problem_data.mu_lame (x{:}), 1, msh_side.nqn, msh_side.nel);
            coeff = (bsxfun (@times, eps_normal, coeff) + p_normal - gside(x{:})).^2;

            if (strcmpi (flag, 'elements'))
              est_level_patch = sum (reshape (sum (coeff, 1), msh_side.nqn, msh_side.nel));
              inds_patch = get_volumetric_indices (iside, msh_patch.nel_dir, elements);
              inds_level = patch_shifting(iptc) + inds_patch;
              [~,~,inds] = intersect (inds_level, hmsh.active{ilev});
              indices = vol_shifting_indices(ilev) + inds;
              est(indices) = est_level_patch;
%            elseif (strcmpi (flag, 'functions'))
%              msh_side = msh_evaluate_element_list (msh_patch.boundary(iside), elements);
%              sp_patch = hspace.space_of_level(ilev).sp_patch{iptc};
%              sp_bnd = sp_evaluate_element_list (sp_patch.boundary(iside), msh_side, 'value', true);
%              est_level = op_f_v (sp_bnd, msh_side, coeff);
%              dofs_bnd = hspace.boundary.space_of_level(ilev).gnum{iptc_bnd};
%              Csub = hspace.boundary.Csub{ilev}(dofs_bnd,:);
%              dofs = hspace.boundary.dofs(1:ndofs(ilev));
%              est(dofs) = est(dofs) + Csub.' * est_level;
            end
          end
        end          
      end
    end
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function est = compute_jump_terms (u, press, hmsh, hspace, hspace_p, mu_lame, flag)
% Compute the jump terms for multipatch geometries, when marking by functions

  if (strcmpi (flag, 'elements'))
    est = zeros (hmsh.nel, 1);
  elseif (strcmpi (flag, 'functions'))
    est = zeros (hspace.ndof + hspace_p.ndof, 1);
  end

  interfaces = hspace.space_of_level(1).interfaces;

  for iref = 1:numel(interfaces)
% Generate an auxiliary hierarchical mesh, such that the elements on the
%  interface coincide from each side. This is required for integration
    [hmsh_aux, interface_elements, interface_active_elements] = hmsh_refined_mesh_for_interface (hmsh, interfaces(iref));
    if (hmsh.nel ~= hmsh_aux.nel)
      hspace_aux = hspace_in_finer_mesh (hspace, hmsh, hmsh_aux);
      hspace_p_aux = hspace_in_finer_mesh (hspace_p, hmsh, hmsh_aux);
    else
      hspace_aux = hspace;
      hspace_p_aux = hspace_p;
    end

% Compute the integral of the jump of the normal derivative at the interface
    if (strcmpi (flag, 'elements'))
      est_edges = integral_term_by_elements (u, press, hmsh_aux, hspace_aux, hspace_p_aux, interfaces(iref), interface_elements, mu_lame);
      est(interface_active_elements(1,:)) = est(interface_active_elements(1,:)) + est_edges;
      est(interface_active_elements(2,:)) = est(interface_active_elements(2,:)) + est_edges;
    elseif (strcmpi (flag, 'functions'))
      est = est + integral_term_by_functions (u, hmsh_aux, hspace_aux, interfaces(iref), interface_elements, lambda_lame, mu_lame);
    end
  end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function est = integral_term_by_functions (u, press, hmsh, hspace, hspace_p, interface, interface_elements, coeff_mu)
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
      ndofp_until_lev = sum (hspace_p.ndof_per_level(1:lev));
      Nelem = cumsum ([0, hmsh.mesh_of_level(lev).nel_per_patch]);
      u_lev = hspace.Csub{lev} * u(1:ndof_until_lev);
      press_lev = hspace_p.Csub{lev} * press(1:ndofp_until_lev);

      grad_dot_normal = cell (2, 1);
      for ii = [2 1] % This ordering allows me to keep the variables from the first (master) side
        gnum = hspace.space_of_level(lev).gnum{patch(ii)};
        gnum_p = hspace_p.space_of_level(lev).gnum{patch(ii)};
        msh_patch_lev = hmsh.mesh_of_level(lev).msh_patch{patch(ii)};
 
% Set of active elements on the patch that are adjacent to the interface
% The numbering in interface_elements is already ordered to make them automatically coincide
        element_list = get_boundary_indices (side(ii), msh_patch_lev.nel_dir, interface_elements{lev}{ii}-Nelem(patch(ii)));

        msh_side_int = msh_boundary_side_from_interior (msh_patch_lev, side(ii));

        msh_side = msh_eval_boundary_side (msh_patch_lev, side(ii), element_list);
        msh_side_aux = msh_evaluate_element_list (msh_side_int, element_list);
 
        sp_bnd = hspace.space_of_level(lev).sp_patch{patch(ii)}.constructor (msh_side_int);
        spu = sp_evaluate_element_list (sp_bnd, msh_side_aux, 'gradient', true);
        sp_bnd = hspace_p.space_of_level(lev).sp_patch{patch(ii)}.constructor (msh_side_int);
        spp = sp_evaluate_element_list (sp_bnd, msh_side_aux, 'value', true);

        grad = sp_eval_msh (u_lev(gnum), spu, msh_side_aux, 'gradient');
        gradt = permute (grad, [2 1 3 4]);
        normal = reshape (msh_side.normal, 1, [], msh_side.nqn, msh_side.nel);
        eps_normal = reshape (sum (bsxfun (@times, grad + gradt, normal), 2), [], msh_side.nqn, msh_side.nel);

        press_val = reshape (sp_eval_msh (press_lev(gnum_p), spp, msh_side_aux, 'value'), 1, msh_side.nqn, msh_side.nel);
        p_normal = reshape (bsxfun (@times, press_val, msh_side.normal), [], msh_side.nqn, msh_side.nel);

        for idim = 1:hmsh.rdim
          x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
        end
        coeffs = reshape (coeff_mu (x{:}), 1, msh_side.nqn, msh_side.nel);
        grad_dot_normal{ii} = bsxfun (@times, eps_normal, coeffs) + p_normal;
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
function est_edges = integral_term_by_elements (u, press, hmsh, hspace, hspace_p, interface, interface_elements, coeff_mu)
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
      ndofp_until_lev = sum (hspace_p.ndof_per_level(1:lev));
      Nelem = cumsum ([0, hmsh.mesh_of_level(lev).nel_per_patch]);
      u_lev = hspace.Csub{lev} * u(1:ndof_until_lev);
      press_lev = hspace_p.Csub{lev} * press(1:ndofp_until_lev);

      grad_dot_normal = cell (2, 1);
      for ii = [2 1] % This ordering allows me to keep the variables from the first (master) side
        gnum = hspace.space_of_level(lev).gnum{patch(ii)};
        gnum_p = hspace_p.space_of_level(lev).gnum{patch(ii)};
        msh_patch_lev = hmsh.mesh_of_level(lev).msh_patch{patch(ii)};
 
% Set of active elements on the patch that are adjacent to the interface
% The numbering in interface_elements is already ordered to make them automatically coincide
        element_list = get_boundary_indices (side(ii), msh_patch_lev.nel_dir, interface_elements{lev}{ii}-Nelem(patch(ii)));

        msh_side_int = msh_boundary_side_from_interior (msh_patch_lev, side(ii));

        msh_side = msh_eval_boundary_side (msh_patch_lev, side(ii), element_list);
        msh_side_aux = msh_evaluate_element_list (msh_side_int, element_list);
 
        sp_bnd = hspace.space_of_level(lev).sp_patch{patch(ii)}.constructor (msh_side_int);
        spu = sp_evaluate_element_list (sp_bnd, msh_side_aux, 'gradient', true);
        sp_bnd = hspace_p.space_of_level(lev).sp_patch{patch(ii)}.constructor (msh_side_int);
        spp = sp_evaluate_element_list (sp_bnd, msh_side_aux, 'value', true);

        grad = sp_eval_msh (u_lev(gnum), spu, msh_side_aux, 'gradient');
        gradt = permute (grad, [2 1 3 4]);
        normal = reshape (msh_side.normal, 1, [], msh_side.nqn, msh_side.nel);
        eps_normal = reshape (sum (bsxfun (@times, grad + gradt, normal), 2), [], msh_side.nqn, msh_side.nel);
        
        press_val = reshape (sp_eval_msh (press_lev(gnum_p), spp, msh_side_aux, 'value'), 1, msh_side.nqn, msh_side.nel);
        p_normal = reshape (bsxfun (@times, press_val, msh_side.normal), [], msh_side.nqn, msh_side.nel);

        for idim = 1:hmsh.rdim
          x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
        end
        coeffs = reshape (coeff_mu (x{:}), 1, msh_side.nqn, msh_side.nel);
        grad_dot_normal{ii} = bsxfun (@times, eps_normal, coeffs) + p_normal;
      end

      % Reorder quadrature points, to consider the relative orientation of the patches
      grad_dot_normal{2} = reorder_quad_points (grad_dot_normal{2}, interface, msh_side.nqn_dir);
      w = msh_side.quad_weights .* msh_side.jacdet;
      dudn_jump = grad_dot_normal{1} + grad_dot_normal{2};
      dudn_jump = reshape (sum (dudn_jump.^2, 1), msh_side.nqn, msh_side.nel);
      
      est_edges = cat (2, est_edges, sum (w.*dudn_jump));
    end
  end
  est_edges = est_edges(:);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function field = reorder_quad_points (field, interface, nqn_dir)
% Reorder quadrature points of adjacent patches, to get a corresponding numbering
% The indices for quadrature points are always in second-last position

  size_of_field = size (field);

  indices = num2cell (repmat (':', 1, numel (size_of_field)));

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

  % This check is necessary when there is one single element (the last dimension is automatically removed)
  if (size_of_field(end-1) == prod(nqn_dir))
    indices{end-1} = qpoints;
  else
    indices{end} = qpoints;
  end
  field = field(indices{:});

end
