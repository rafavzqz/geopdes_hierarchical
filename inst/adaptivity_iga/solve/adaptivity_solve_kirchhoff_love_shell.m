% ADAPTIVITY_SOLVE_KIRCHHOFF_LOVE_SHELL: Solve the Kirchhoff-Love shell problem with hierarchical splines.
%
% USAGE:
%
% u = adaptivity_solve_kirchhoff_love_shell (hmsh, hspace, problem_data)
%
% INPUT:
%
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the vector-valued space of hierarchical splines (see hierarchical_space)
%   problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - drchlt_sides: sidues with Dirichlet boundary condition
%    - E_coeff:      Young's modulus
%    - nu_coeff:     Poisson ratio
%    - thickness:    thickness of the shell
%    - f:            source term
%
% OUTPUT:
%
%   u: computed degrees of freedom
%
% Copyright (C) 2023 Rafael Vazquez
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

function [u, energy] = adaptivity_solve_kirchhoff_love_shell (hmsh, hspace, problem_data)

Kmat = op_KL_shells_hier (hspace, hspace, hmsh, problem_data.E_coeff, problem_data.nu_coeff, problem_data.thickness);
rhs = op_f_v_hier (hspace, hmsh, problem_data.f);

% Apply boundary conditions
% TODO : so far, only homogeneous boundary conditions in every component are implemented
u = zeros (hspace.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (hspace, hmsh, @(x,y,z,ind) zeros([3, size(x)]), problem_data.drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:hspace.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - Kmat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u(int_dofs) = Kmat(int_dofs, int_dofs) \ rhs(int_dofs);
energy = 0.5 * u(int_dofs).' * Kmat(int_dofs, int_dofs) * u(int_dofs);

end