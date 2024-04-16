% ADAPTIVITY_SOLVE_KIRCHHOFF_LOVE_SHELL_MP_C1: solve the Kirchhoff-Love shell problem 
%  for hierarhical multipatch C^1 splines.
%
% USAGE:
%
%  [u, energy] = adaptivity_solve_kirchhoff_love_shell_mp_C1 (hmsh, hspace, problem_data)
%
% INPUT:
%
%  hmsh:         object representing the hierarchical mesh (see hierarchical_mesh_mp)
%  hspace:       object representing the space of hierarchical splines (see hierarchical_space_mp_C1)
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - drchlt_sides: sides with Dirichlet boundary condition (for now, only homogeneous conditions in every direction)
%    - components:   vector components on which to apply the zero displacement condition, a cell-array with the length of drchlt_sides
%    - E_coeff:      function handle for Young's modulus
%    - nu_coeff:     function handle for Poisson's ratio
%    - thickness:    scalar value, thickness of the shell
%    - f:            source term, distributed load
%
% OUTPUT:
%
%  u:        the computed degrees of freedom
%  energy:   energy of the system (sqrt(u' * K * u)), for postprocessing.
%
% Copyright (C) 2022, 2023 Cesare Bracco, Andrea Farahat, Rafael Vazquez
% Copyright (C) 2024 Rafael Vazquez
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

function [u, energy] = adaptivity_solve_kirchhoff_love_shell_mp_C1 (hmsh, hspace, problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end

Kmat = op_KL_shells_hier (hspace, hspace, hmsh, E_coeff, nu_coeff, thickness);
rhs = op_f_v_hier_vector (hspace, hmsh, f);
% Apply zero rotation with Nitsche method
if (exist ('rotation_sides', 'var') && ~isempty(rotation_sides))
  Kmat = Kmat + sp_nitsche_KL_rotation (hspace, hmsh, rotation_sides, E_coeff, nu_coeff, thickness, method_data.penalty_coeff);
end

% Apply boundary conditions
% TODO : so far, only homogeneous boundary conditions in every component are implemented
[~, drchlt_dofs, kernel_dofs] = sp_drchlt_C1_shells (hspace, hmsh, drchlt_sides);

ndof = hmsh.rdim * hspace.ndof;
u = zeros (ndof, 1);
u(drchlt_dofs) = 0; % TODO: use u_drchlt;

int_dofs = setdiff (1:ndof, drchlt_dofs);
add_dofs = kernel_dofs.quasi_interior_dofs; %this will contain the "boundary" vertex dofs which have been removed from drchlt_dofs
add_dofs = [add_dofs, hspace.ndof+add_dofs, 2*hspace.ndof+add_dofs];

% We assemble the (pieces of the) stiffness matrix, the rhs (and its correction taking 
% into account the Dirichlet conditions), and the basis change matrix (we will need it 
% to go from the basis with kernel vectors obtained when examining the Dirichlet conditions 
% to the usual basis)
vertex_dofs = kernel_dofs.all_vertex_dofs;
vertex_dofs = [vertex_dofs, hspace.ndof+vertex_dofs, 2*hspace.ndof+vertex_dofs];

B_change_vector = blkdiag(kernel_dofs.B_change_local, kernel_dofs.B_change_local, kernel_dofs.B_change_local);
Kmat(:,add_dofs) = Kmat(:,vertex_dofs) * B_change_vector;
Kmat(add_dofs,:) = B_change_vector.' * Kmat(vertex_dofs,:);
rhs(add_dofs) = B_change_vector.' * rhs(vertex_dofs);

% Solve the linear system
u(int_dofs) = Kmat(int_dofs, int_dofs) \ rhs(int_dofs);
% Compute the energy (only valid for homogeneous b.c.)
energy = 0.5 * u(int_dofs).' * Kmat(int_dofs, int_dofs) * u(int_dofs);

% Switching to the usual basis using the local matrix for the vertex dofs
u_old = u(setdiff(vertex_dofs, add_dofs)); % Coefficients of the vertex functions that were already in the old basis
u(vertex_dofs) = B_change_vector * u(add_dofs);
u(setdiff(vertex_dofs, add_dofs)) = u(setdiff(vertex_dofs, add_dofs)) + u_old;

end
