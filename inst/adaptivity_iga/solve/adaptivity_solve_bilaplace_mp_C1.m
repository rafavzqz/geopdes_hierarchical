% ADAPTIVITY_SOLVE_BILAPLACE_MPC1: Solve a 2d bilaplace problem with a
% variational formulation using the laplacian
%
% The function solves the bilaplacian problem
%
%      epsilon laplace(laplace(u)) = f    in Omega = F((0,1)^2)
%                            du/dn = 0    on Gamma_D (Rotation)
%                                u = 0    on Gamma_D (Transverse displacement)
%           m = epsilon d^2 u/dn^2 = 0    on Gamma_N (Distributed double forces)
%
% with a variational formulation using the Hessian matrix, that is
%
%       epsilon ( laplace(u), laplace(v) ) = (f,v)
%
% USAGE:
%
%  [geometry, msh, space, u] = solve_bilaplace_mpC1 (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - drchlt_sides: sides with Dirichlet boundary condition (always homogeneous)
%    - c_diff:       constant physical parameter (epsilon in the equation)
%    - f:            source term
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space:    space object that defines the discrete space (see sp_scalar)
%  u:        the computed degrees of freedom
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, 2013 Rafael Vazquez
% Copyright (C) 2013, Marco Pingaro
% Copyright (C) 2020-2023, Cesare Bracco, Andrea Farahat, Rafael Vazquez
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

function u = adaptivity_solve_bilaplace_mp_C1 (hmsh, hspace, problem_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end

% Assemble the matrices
stiff_mat = op_laplaceu_laplacev_hier (hspace, hspace, hmsh, c_diff);  %use hierarchical functions here
rhs       = op_f_v_hier (hspace, hmsh, f);

% Apply essential conditions for the value and the normal derivative
%  The kernel_dofs are used to avoid linear dependencies on the boundary
u = zeros (hspace.ndof, 1);
if (isfield(problem_data, 'graduex') && isfield(problem_data, 'uex'))
    [u_drchlt, drchlt_dofs, kernel_dofs] = sp_bilaplacian_drchlt_C1_exact (hspace, hmsh, drchlt_sides, uex, graduex);
else
    [u_drchlt, drchlt_dofs, kernel_dofs] = sp_bilaplacian_drchlt_C1 (hspace, hmsh, drchlt_sides, h, g);
end
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:hspace.ndof, drchlt_dofs);
add_dofs = kernel_dofs.quasi_interior_dofs; % This will contain the "boundary" vertex dofs which have been removed from drchlt_dofs

% We assemble the stiffness matrix and the rhs taking into account the kernel dofs
vertex_dofs = kernel_dofs.all_vertex_dofs;

stiff_mat(:,add_dofs) = stiff_mat(:,vertex_dofs) * kernel_dofs.B_change_local;
stiff_mat(add_dofs,:) = kernel_dofs.B_change_local.' * stiff_mat(vertex_dofs,:);
rhs(add_dofs) = kernel_dofs.B_change_local.' * rhs(vertex_dofs);

rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);

% Switching to the usual basis using the local matrix for the vertex dofs
u_old = u(setdiff(kernel_dofs.all_vertex_dofs, add_dofs)); % Coefficients of the vertex functions that were already in the old basis
u(kernel_dofs.all_vertex_dofs) = kernel_dofs.B_change_local * u(kernel_dofs.quasi_interior_dofs);
u(setdiff(kernel_dofs.all_vertex_dofs, add_dofs)) = u(setdiff(kernel_dofs.all_vertex_dofs, add_dofs)) + u_old;

end
