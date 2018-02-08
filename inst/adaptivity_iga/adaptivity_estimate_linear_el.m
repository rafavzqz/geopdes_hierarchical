% ADAPTIVITY_ESTIMATE_LINEAR_EL: Computation of a posteriori error indicators for linear elasticity problem, using globally smooth (C^1) hierarchical spaces.
%
% We consider the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^n)
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% USAGE:
%
%   est = adaptivity_estimate_lINEAR_EL (u, hmsh, hspace, problem_data, adaptivity_data)
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


function est = adaptivity_estimate_linear_el (u, hmsh, hspace, problem_data, adaptivity_data)

if (isfield(adaptivity_data, 'C0_est'))
    C0_est = adaptivity_data.C0_est;
else
    C0_est = 1;
end

[ders2, F] = hspace_eval_hmsh (u, hspace, hmsh, 'hessian');

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
for h = 1:hspace.ncomp
    partials_a = 0;
    partials_b = 0;
    for t = 1:hspace.ncomp
        if (t ~= h)
            partials_a = partials_a + reshape (ders2(t,h,t,:,:), [], hmsh.nel); %mixed derivatives of all components (except h-th component)
            partials_b = partials_b + reshape (ders2(h,t,t,:,:), [], hmsh.nel); %second derivatives of h-th component (except w.r.t. h-th variable)
        end
    end
    divergence(h,:,:) = (2*mu+lambda)*reshape(ders2(h,h,h,:,:), [], hmsh.nel) + (mu+lambda)*partials_a + mu*partials_b;
end
aux = (valf + divergence).^2;  %residual
aux = reshape (sum(aux), [], hmsh.nel);

switch adaptivity_data.flag
    case 'elements',
        w = [];
        h = [];
        for ilev = 1:hmsh.nlevels
            if (hmsh.msh_lev{ilev}.nel ~= 0)
                w = cat (2, w, hmsh.msh_lev{ilev}.quad_weights .* hmsh.msh_lev{ilev}.jacdet);
                h = cat (1, h, hmsh.msh_lev{ilev}.element_size(:));
            end
        end
        h = h * sqrt (hmsh.ndim);
        
        est = sqrt (sum (aux.*w));
        est = C0_est*h.*est(:);

    % Jump terms, only computed for multipatch geometries
        if (isa (hmsh, 'hierarchical_mesh_mp') && hmsh.npatch > 1)
          jump_est = compute_jump_terms (u, hmsh, hspace, problem_data.lambda_lame, problem_data.mu_lame, adaptivity_data.flag);
          est = est + h .* jump_est;
        end
        est = C0_est * sqrt (est);
        
    case 'functions',
        ms = zeros (hmsh.nlevels, 1);
        for ilev = 1:hmsh.nlevels
            if (hmsh.msh_lev{ilev}.nel ~= 0)
                ms(ilev) = max (hmsh.msh_lev{ilev}.element_size);
            else
                ms(ilev) = 0;
            end
        end
        ms = ms * sqrt (hmsh.ndim);
        
        Nf = cumsum ([0; hspace.ndof_per_level(:)]);
        dof_level = zeros (hspace.ndof, 1);
        for lev = 1:hspace.nlevels
            dof_level(Nf(lev)+1:Nf(lev+1)) = lev;
        end
        coef = ms(dof_level).*sqrt(hspace.coeff_pou(:));
        
        est = zeros(hspace.ndof,1);
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
        est = C0_est * coef .* sqrt(est);
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function est = compute_jump_terms (u, hmsh, hspace, lambda_lame, mu_lame, flag)
% Compute the jump terms for multipatch geometries, when marking by functions

  if (strcmpi (flag, 'elements'))
    est = zeros (hmsh.nel, 1);
%   elseif (strcmpi (flag, 'functions'))
%     est = zeros (hspace.ndof, 1);
  end

  interfaces = hspace.space_of_level(1).interfaces;

  for iref = 1:numel(interfaces)
% Generate an auxiliary hierarchical mesh, such that the elements on the
%  interface coincide from each side. This is required for integration
    [hmsh_aux, interface_elements, interface_active_elements] = generate_auxiliary_mesh (hmsh, interfaces(iref));
    if (hmsh.nel ~= hmsh_aux.nel)
      hspace_aux = hspace_in_finer_mesh (hspace, hmsh, hmsh_aux);
    else
      hspace_aux = hspace;
    end

% Compute the integral of the jump of the normal derivative at the interface
    if (strcmpi (flag, 'elements'))
      est_edges = integral_term_by_elements (u, hmsh_aux, hspace_aux, interfaces(iref), interface_elements, lambda_lame, mu_lame);
      est(interface_active_elements(1,:)) = est(interface_active_elements(1,:)) + est_edges;
      est(interface_active_elements(2,:)) = est(interface_active_elements(2,:)) + est_edges;
    elseif (strcmpi (flag, 'functions'))
      est = est + integral_term_by_functions (u, hmsh_aux, hspace_aux, interfaces(iref), interface_elements, lambda_lame, mu_lame);
    end
  end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hmsh_aux, interface_elements, adjacent_elements_to_edge] = generate_auxiliary_mesh (hmsh, interface)
% Generate an auxiliary (refined) hierarchical mesh, such that two adjacent elements
%  on the interface are active on both patches
%
% OUTPUT
%    hmsh_aux:           refined hierarchical mesh
%    interface_elements: for each level, and for each side of the interface,
%           indices of the adjacent active elements in hmsh_aux
%    adjacent_elements_to_edge: array of size (2, nedges). For each edge of the interface,
%           and for each side, global indices of the adjacent active elements (in hmsh)

  patch(1) = interface.patch1;
  patch(2) = interface.patch2;
  side(1) = interface.side1;
  side(2) = interface.side2;
  
  iedge = 0;
  hmsh_aux = hmsh;
  interface_elements = cell (hmsh.nlevels, 1);
  
  for lev = 1:hmsh.nlevels
    marked = cell (hmsh.nlevels, 1);
    Nelem = cumsum ([0, hmsh_aux.mesh_of_level(lev).nel_per_patch]);
    for ii = 1:2
      msh_patch_lev = hmsh_aux.mesh_of_level(lev).msh_patch{patch(ii)};
      nel_dir = msh_patch_lev.nel_dir;
%    ind = [1 1 2 2 3 3] in 3D, ind = [1 1 2 2] in 2D
%    ind2 = [2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind = [2 2 1 1] in 2D;
      ind = ceil (side(ii)/2);
      ind2 = setdiff (1:hmsh.ndim, ind);
      subindices = arrayfun (@(x) 1:x, nel_dir, 'UniformOutput', false);
      if (mod (side(ii), 2) == 1)
        subindices{ind} = 1;
      else
        subindices{ind} = nel_dir(ind);
      end
      [subindices{:}] = ndgrid (subindices{:});
      elems{ii} = reshape (sub2ind ([nel_dir, 1], subindices{:}), [nel_dir(ind2), 1]);
    end
    elems{1} = elems{1}(:)';
    elems{2} = reorder_elements (elems{2}, interface, nel_dir(ind2));
    [active_elements{1}, pos1] = ismember (elems{1}+Nelem(patch(1)), hmsh_aux.active{lev});
    [active_elements{2}, pos2] = ismember (elems{2}+Nelem(patch(2)), hmsh_aux.active{lev});
    
    interface_indices = active_elements{1} & active_elements{2};
    interface_elements{lev}{1} = hmsh_aux.active{lev}(pos1(interface_indices));
    interface_elements{lev}{2} = hmsh_aux.active{lev}(pos2(interface_indices));

    indices1 = (active_elements{1} & ~active_elements{2});
    indices2 = (active_elements{2} & ~active_elements{1});
    indices = union (pos1(indices1), pos2(indices2));
    marked{lev} = hmsh_aux.active{lev}(indices);
    hmsh_aux = hmsh_refine (hmsh_aux, marked);

    iedge_aux = iedge;
    Nelem_level = cumsum ([0 hmsh.nel_per_level]);
    for ii = 1:2
      iedge = iedge_aux;
      int_elems = interface_elements{lev}{ii};
      for iel = 1:numel(int_elems)
        elem_lev = int_elems(iel);
        iedge = iedge + 1;
        flag = 0;
        levj = lev;
        while (~flag)
          [is_active, pos] = ismember (elem_lev, hmsh.active{levj});
          if (is_active)
            adjacent_elements_to_edge(ii,iedge) = Nelem_level(levj) + pos;
            flag = 1;
          else
            elem_lev = hmsh_get_parent (hmsh, levj, elem_lev);
            flag = 0;
            levj = levj-1;
          end
        end
      end
    end

  end  

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function est_edges = integral_term_by_elements (u, hmsh, hspace, interface, interface_elements, coeff_lambda, coeff_mu)
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
        spp = sp_evaluate_element_list (sp_bnd, msh_side_aux, 'gradient', true, 'divergence', true);

        [grad_and_div] = sp_eval_msh (u_lev(gnum), spp, msh_side_aux, {'gradient', 'divergence'});
        [grad, div] = deal (grad_and_div{:});
        gradt = permute (grad, [2 1 3 4]);
        normal = reshape (msh_side.normal, 1, [], msh_side.nqn, msh_side.nel);
        eps_normal = reshape (sum (bsxfun (@times, grad + gradt, normal), 2), [], msh_side.nqn, msh_side.nel);
        
        div = reshape (div, 1, msh_side.nqn, msh_side.nel);
        div_normal = reshape (bsxfun (@times, div, msh_side.normal), [], msh_side.nqn, msh_side.nel);

        for idim = 1:hmsh.rdim
          x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
        end
        coeffs = reshape (coeff_mu (x{:}), 1, msh_side.nqn, msh_side.nel);
        coeffs_lambda = reshape (coeff_lambda (x{:}), 1, msh_side.nqn, msh_side.nel);
        grad_dot_normal{ii} = bsxfun (@times, eps_normal, coeffs) + bsxfun (@times, div_normal, coeffs_lambda);
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
function elem = reorder_elements (elem, interface, nel_dir)
% Reorder elements of adjacent patches, to get a corresponding numbering
  ndim = numel (nel_dir) + 1;
  if (ndim == 2)
    if (interface.ornt == -1)
      elem = fliplr (elem(:)');
    else
      elem = elem(:)';
    end
  elseif (ndim == 3)
    if (interface.flag == -1)
      elem = elem';
    end
    if (interface.ornt1 == -1)
      elem = flipud (elem);
    end
    if (interface.ornt2 == -1)
      elem = fliplr (elem);
    end
    elem = elem(:)';
  end
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
