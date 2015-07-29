% ADAPTIVITY_SOLVE: assemble and solve the linear system for Laplacian problem, using hierarchical spaces.
%
% The function solves the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^n)
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% USAGE:
%
% u = adaptivity_solve (hmsh, hspace, method_data)
%
% INPUT:
%
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the space of hierarchical splines (see hierarchical_space)
%   problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_diff:       diffusion coefficient (epsilon in the equation)
%    - f:            source term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
%
% OUTPUT:
%
%   u: computed degrees of freedom
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

function u = adaptivity_solve (hmsh, hspace, problem_data)

stiff_mat = op_gradu_gradv_hier (hspace, hspace, hmsh, problem_data.c_diff);
rhs = op_f_v_hier (hspace, hmsh, problem_data.f);
% mass = op_u_v_hier (hspace, hspace, hmsh, problem_data.c_diff);

% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX TO BE TESTED
% Apply Neumann boundary conditions 
for iside = problem_data.nmnn_sides
  if (hmsh.ndim > 1)
% Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
    gside = @(varargin) problem_data.g(varargin{:},iside);
    dofs = hspace.boundary(iside).dofs;
    rhs(dofs) = rhs(dofs) + op_f_v_hier (hspace.boundary(iside), hmsh.boundary(iside), gside);
  else
    if (iside == 1)
      x = hmsh.mesh_of_level(1).breaks{1}(1);
    else
      x = hmsh.mesh_of_level(1).breaks{1}(end);
    end
    sp_side = hspace.boundary(iside);
    rhs(sp_side.dofs) = rhs(sp_side.dofs) + problem_data.g(x,iside);
  end
end

% Apply Dirichlet boundary conditions
u = zeros (hspace.ndof, 1);
[u_dirichlet, dirichlet_dofs] = hsp_drchlt_l2_proj (hspace, hmsh, problem_data.h, problem_data.drchlt_sides);
u(dirichlet_dofs) = u_dirichlet;

int_dofs = setdiff (1:hspace.ndof, dirichlet_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, dirichlet_dofs)*u(dirichlet_dofs);

% Solve the linear system
u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);

fprintf('Total DOFs: %d - Used DOFs: %d (%f seconds)\n', hspace.ndof, numel (int_dofs));
