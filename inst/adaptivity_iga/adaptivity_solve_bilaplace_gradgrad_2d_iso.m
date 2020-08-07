% ADAPTIVITY_SOLVE_BILAPLACE_2D_GRADGRAD_ISO: Solve a 2d bilaplace problem
% with a variational formulation using the gradient of the gradient with
% hierarchical splines
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
%       epsilon ( grad(grad(u)), grad(grad(v)) ) = (f,v)
%
% where grad(grad(u)) = d^2 u/ dx_i dx_j
%
% USAGE:
%
%  [geometry, hmsh, hspace, u] = adaptivity_solve_bilaplace_2d_gradgrad_iso (problem_data, method_data)
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
% See also EX_KIRCHHOFF_RECTANGULAR_PLATE, EX_KIRCHHOFF_RECTANGULAR_PLATE_CLAMPED 
% EX_KIRCHHOFF_CIRCULAR_PLATE, EX_KIRCHHOFF_CIRCULAR_PLATE_CLAMPED for some examples.
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, 2013 Rafael Vazquez
% Copyright (C) 2013, Marco Pingaro
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

function u = adaptivity_solve_bilaplace_gradgrad_2d_iso (hspace, hmsh, problem_data)

% Assemble the matrices
stiff_mat = op_gradgradu_gradgradv_hier (hspace, hspace, hmsh, problem_data.D);
rhs       = op_f_v_hier (hspace, hmsh, problem_data.f);

if (isfield(problem_data, 'point_load'))
    contribution_to_rhs = op_point_load_hier(hspace, hmsh, problem_data.point_load, problem_data.local_point);
    rhs = rhs + contribution_to_rhs;
end

% Apply simply supported conditions
u = zeros (hspace.ndof, 1);
drchlt_dofs_u = []; drchlt_dofs_r = []; dirichlet_dofs = [];
for iside = problem_data.simply_supported_sides
  drchlt_dofs_u = union (drchlt_dofs_u, hspace.boundary(iside).dofs);
end

% Apply Neumann boundary conditions
for iside = problem_data.nmnn_sides
% Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
  gside = @(varargin) problem_data.g(varargin{:},iside);
  dofs = hspace.boundary(iside).dofs;
  rhs(dofs) = rhs(dofs) + op_f_v_hier (hspace.boundary(iside), hmsh.boundary(iside), gside);
end

% Apply inhomogeneous weak conditions on moments
for iside = problem_data.prescribed_moment_sides
  N_rhs = op_gradv_n_f_hier (hspace, hmsh, problem_data.M, iside);
  dofs_moment = hspace.boundary(iside).dofs;
  adj_dofs = hspace.boundary(iside).adjacent_dofs;
  rhs(dofs_moment) = rhs(dofs_moment) + N_rhs(dofs_moment);
  rhs(adj_dofs) = rhs(adj_dofs) + N_rhs(adj_dofs);
end

% Apply pressure conditions
for iside = problem_data.press_sides
  pside = @(varargin) problem_data.p(varargin{:},iside);
  dofs = hspace.boundary(iside).dofs;
  rhs(dofs) = rhs(dofs) - op_pn_v_hier (hspace.boundary(iside), hmsh.boundary(iside), pside);
end

% % Apply Dirichlet boundary conditions
[u_dirichlet, dirichlet_dofs] = sp_drchlt_l2_proj (hspace, hmsh, problem_data.h, problem_data.drchlt_sides);
u(dirichlet_dofs) = u_dirichlet;

% Apply clamped conditions
for iside = problem_data.clamped_sides
  drchlt_dofs_u = union (drchlt_dofs_u, hspace.boundary(iside).dofs);
  drchlt_dofs_r = union (drchlt_dofs_r, hspace.boundary(iside).adjacent_dofs);
end
dirichlet_dofs_clamped = union (drchlt_dofs_u, drchlt_dofs_r);

dirichlet_dofs = union (dirichlet_dofs, dirichlet_dofs_clamped);

int_dofs = setdiff (1:hspace.ndof, dirichlet_dofs);

rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, dirichlet_dofs)*u(dirichlet_dofs);

% Solve the linear system
u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);
end
