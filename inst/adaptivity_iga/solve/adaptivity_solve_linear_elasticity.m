% ADAPTIVITY_SOLVE_LINEAR_ELASTICITY: Solve a linear elasticity problem with hierarchical splines.
%
% The function solves the linear elasticity problem
%
%      - div (sigma(u)) = f    in Omega = F((0,1)^n)
%      sigma(u) \cdot n = g    on Gamma_N
%                     u = h    on Gamma_D
%
% with   sigma(u) = mu*(grad(u) + grad(u)^t) + lambda*div(u)*I.
%
%   u:          displacement vector
%   sigma:      Cauchy stress tensor
%   lambda, mu: Lame' parameters
%   I:          identity tensor
%
% USAGE:
%
% u = adaptivity_solve_linear_elasticity (hmsh, hspace, problem_data)
%
% INPUT:
%
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the vector-valued space of hierarchical splines (see hierarchical_space)
%   problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - press_sides:  sides with pressure boundary condition (may be empty)
%    - symm_sides:   sides with symmetry boundary condition (may be empty)
%    - c_diff:       diffusion coefficient (epsilon in the equation)
%    - lambda_lame:  first Lame' parameter
%    - mu_lame:      second Lame' parameter
%    - f:            source term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
%
% OUTPUT:
%
%   u: computed degrees of freedom
%
% Copyright (C) 2015 Eduardo M. Garau
% Copyright (C) 2015-2024 Rafael Vazquez
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

function u = adaptivity_solve_linear_elasticity (hmsh, hspace, problem_data)

mat    = op_su_ev_hier (hspace, hspace, hmsh, problem_data.lambda_lame, problem_data.mu_lame); 
rhs    = op_f_v_hier (hspace, hmsh, problem_data.f);

% Apply Neumann boundary conditions
if (~isfield (struct (hmsh), 'npatch')) % Single patch case
  for iside = problem_data.nmnn_sides
% Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
    gside = @(varargin) problem_data.g(varargin{:},iside);
    dofs = hspace.boundary(iside).dofs;
    rhs(dofs) = rhs(dofs) + op_f_v_hier (hspace.boundary(iside), hmsh.boundary(iside), gside);
  end
else % Multipatch case
  boundaries = hmsh.mesh_of_level(1).boundaries;
  Nbnd = cumsum ([0, boundaries.nsides]);
  for iref = problem_data.nmnn_sides
    iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    gref = @(varargin) problem_data.g(varargin{:},iref);
    rhs_nmnn = op_f_v_hier (hspace.boundary, hmsh.boundary, gref, iref_patch_list);
    rhs(hspace.boundary.dofs) = rhs(hspace.boundary.dofs) + rhs_nmnn;
  end
end

% Apply pressure conditions
if (~isfield (struct (hmsh), 'npatch')) % Single patch case
  for iside = problem_data.press_sides
    error ('Press sides are not implemented yet')
    pside = @(varargin) problem_data.p(varargin{:},iside);
    dofs = hspace.boundary(iside).dofs;
    rhs(dofs) = rhs(dofs) - op_pn_v_hier (hspace.boundary(iside), hmsh.boundary(iside), pside);
  end
else % Multipatch case
  boundaries = hmsh.mesh_of_level(1).boundaries;
  Nbnd = cumsum ([0, boundaries.nsides]);
  for iref = problem_data.press_sides
    error ('Press sides are not implemented yet')
    iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    pref = @(varargin) problem_data.p(varargin{:},iref);
    rhs_press = op_pn_v_hier (hspace.boundary, hmsh.boundary, pref, iref_patch_list);
    rhs(hspace.boundary.dofs) = rhs(hspace.boundary.dofs) - rhs_press;
  end
end

% Apply symmetry conditions
u = zeros (hspace.ndof, 1);
symm_dofs = [];
if (~isfield (struct (hmsh), 'npatch')) % Single patch case
  for iside = problem_data.symm_sides
    msh_side = hmsh_eval_boundary_side (hmsh, iside);
%     normal_comp = zeros (msh_side.rdim, msh_side.nqn * msh_side.nel);
%     for idim = 1:msh_side.rdim
%       normal_comp(idim,:) = reshape (msh_side.normal(idim,:,:), 1, msh_side.nqn*msh_side.nel);
%     end
    normal_comp = reshape (msh_side.normal, msh_side.rdim, msh_side.nqn*msh_side.nel);

    parallel_to_axes = false;
    for ind = 1:msh_side.rdim
      ind2 = setdiff (1:msh_side.rdim, ind);
      if (all (all (abs (normal_comp(ind2,:)) < 1e-10)))
        symm_dofs = union (symm_dofs, hspace.boundary(iside).dofs(hspace.boundary(iside).comp_dofs{ind}));
        parallel_to_axes = true;
        break
      end
    end
    if (~parallel_to_axes)
      error ('adaptivity_solve_linear_elasticity: We have only implemented the symmetry condition for boundaries parallel to the axes')
    end
  end
else % Multipatch case
  boundaries = hmsh.mesh_of_level(1).boundaries;
  Nbnd = cumsum ([0, boundaries.nsides]);
  symm_dofs = [];
  for iref = problem_data.symm_sides
    iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    for ii = 1:numel(iref_patch_list)
      iptc_bnd = iref_patch_list(ii);
      iptc = boundaries(iref).patches(ii);
      iside = boundaries(iref).faces(ii);
      msh_patch = hmsh.mesh_of_level(1).msh_patch{iptc};
      msh_side = msh_eval_boundary_side (msh_patch, iside);
      normal_comp = reshape (msh_side.normal, msh_side.rdim, msh_side.nqn*msh_side.nel);

      parallel_to_axes = false;
      for ind = 1:msh_side.rdim
        ind2 = setdiff (1:msh_side.rdim, ind);
        if (all (all (abs (normal_comp(ind2,:)) < 1e-10)))
          parallel_to_axes = true;
          break
        end
      end
      if (~parallel_to_axes)
        error ('adaptivity_solve_linear_elasticity: We have only implemented the symmetry condition for boundaries parallel to the axes')
      else
        shifting_indices = cumsum ([0 hspace.boundary.ndof_per_level]);
        symm_dofs_iref = [];
        for ilev = 1:hmsh.boundary.nlevels
          sp_lev = hspace.boundary.space_of_level(ilev);
          global_comp_dofs = sp_lev.gnum{iptc_bnd}(sp_lev.sp_patch{iptc_bnd}.comp_dofs{ind});
          [~,~,bnd_indices] = intersect (global_comp_dofs, hspace.boundary.active{ilev});
          symm_dofs_iref = union (symm_dofs_iref, shifting_indices(ilev) + bnd_indices);
        end
        symm_dofs = union (symm_dofs, hspace.boundary.dofs(symm_dofs_iref));
      end
    end
  end
end  

% Apply Dirichlet boundary conditions
[u_dirichlet, dirichlet_dofs] = sp_drchlt_l2_proj (hspace, hmsh, problem_data.h, problem_data.drchlt_sides);
u(dirichlet_dofs) = u_dirichlet;

int_dofs = setdiff (1:hspace.ndof, union (dirichlet_dofs, symm_dofs));
rhs(int_dofs) = rhs(int_dofs) - mat(int_dofs, dirichlet_dofs)*u_dirichlet;

% Solve the linear system
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

end
