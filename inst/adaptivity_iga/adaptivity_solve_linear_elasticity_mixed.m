% ADAPTIVITY_SOLVE_LINEAR_ELASTICITY_MIXED: Solve a linear elasticity problem with hierarchical splines, using a mixed formulation.
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
% [u, press] = adaptivity_solve_linear_elasticity_mixed (hmsh, hspace, hspace_mul, method_data)
%
% INPUT:
%
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the vector-valued space of hierarchical splines for the displacement (see hierarchical_space)
%   hspace_mul: object representing the scalar-valued space of hierarchical splines for the Lagrange multiplier (see hierarchical_space)
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
% Copyright (C) 2015, 2017 Rafael Vazquez
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

function [vel, press] = adaptivity_solve_linear_elasticity_mixed (hmsh, hspace, hspace_mul, problem_data)

invlam = @(varargin) 1./problem_data.lambda_lame(varargin{:});
A   = op_eu_ev_hier (hspace, hspace, hmsh, problem_data.mu_lame);
B   = op_div_v_q_hier (hspace, hspace_mul, hmsh);
M   = op_u_v_hier (hspace_mul, hspace_mul, hmsh, invlam);
rhs = op_f_v_hier (hspace, hmsh, problem_data.f);

% Apply Neumann boundary conditions
for iside = problem_data.nmnn_sides
% Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
  gside = @(varargin) problem_data.g(varargin{:},iside);
  dofs = hspace.boundary(iside).dofs;
  rhs(dofs) = rhs(dofs) + op_f_v_hier (hspace.boundary(iside), hmsh.boundary(iside), gside);
end

% Apply pressure conditions
for iside = problem_data.press_sides
  pside = @(varargin) problem_data.p(varargin{:},iside);
  dofs = hspace.boundary(iside).dofs;
  rhs(dofs) = rhs(dofs) - op_pn_v_hier (hspace.boundary(iside), hmsh.boundary(iside), pside);
end

% Apply symmetry conditions
vel   = zeros (hspace.ndof, 1);
press = zeros (hspace_mul.ndof, 1);
symm_dofs = [];
for iside = problem_data.symm_sides
  msh_side = hmsh_eval_boundary_side (hmsh, iside);
  normal_comp = zeros (msh_side.rdim, msh_side.nqn * msh_side.nel);
  for idim = 1:msh_side.rdim
    normal_comp(idim,:) = reshape (msh_side.normal(idim,:,:), 1, msh_side.nqn*msh_side.nel);
  end

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
    error ('adaptivity_solve_linear_elasticity_mixed: We have only implemented the symmetry condition for boundaries parallel to the axes')
  end
end

% Apply Dirichlet boundary conditions
[u_dirichlet, dirichlet_dofs] = sp_drchlt_l2_proj (hspace, hmsh, problem_data.h, problem_data.drchlt_sides);
vel(dirichlet_dofs) = u_dirichlet;

int_dofs = setdiff (1:hspace.ndof, union (dirichlet_dofs, symm_dofs));

mat = [ A(int_dofs, int_dofs),  B(:,int_dofs).';
        B(:,int_dofs),          -M];

rhs(int_dofs) = rhs(int_dofs) - A(int_dofs, dirichlet_dofs)*u_dirichlet;
rhs_tot = [rhs(int_dofs); B(:, dirichlet_dofs)*u_dirichlet];

rhs(int_dofs) = rhs(int_dofs) - mat(int_dofs, dirichlet_dofs)*u_dirichlet;

% Solve the linear system
u = mat \ rhs_tot;

press = u(numel(int_dofs)+1:end);
vel(int_dofs) = u(1:numel(int_dofs));
