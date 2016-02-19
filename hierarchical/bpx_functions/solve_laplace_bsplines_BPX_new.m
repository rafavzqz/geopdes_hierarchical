% SOLVE_LAPLACE_2D: Solve a 2d Laplace problem with a T-spline discretization (non-isoparametric approach). 
%
% The function solves the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^2)
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% USAGE:
%
%  [geometry, msh, space, u] = solve_laplace_2d (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_diff:       diffusion coefficient (epsilon in the equation)
%    - f:            source term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
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
%  msh:      mesh object that defines the quadrature rule (see msh_2d)
%  space:    space object that defines the discrete space (see sp_bspline_2d)
%  u:        the computed degrees of freedom
%
% See also EX_LAPLACE_SQUARE for an example.
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, Rafael Vazquez
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

% BPX NOTATION
% Ai:       stiffness matrix of level
% rhs:      right-hand side of level
% int_dofs: internal dofs (usually int_dofs)
% ndof:     total number of dofs (usually ndof)
% new_dofs: basis functions (in int_dofs) that were not present in the previous level
% Pi:       projector between two consecutive levels (from l-1 to l)
% Qi:       projector between the current level and the finest one (from l to nlevels)

function [geometry, msh_1, space_1, u] = ...
              solve_laplace_bsplines_BPX_new (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Construct geometry structure
geometry  = geo_load (geo_name);

[knots, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);

% bpx(1:totlevels) = struct('ndof',0,'dofs',[],'Ai',[],'Qi',[],'Pi',[]);

rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (knots, rule);
msh      = msh_cartesian (knots, qn, qw, geometry);
space    = sp_bspline (knots, degree, msh);

% Assemble the matrices
Ai = op_gradu_gradv_tp (space, space, msh);
rhs       = op_f_v_tp (space, msh, f);

% Apply Dirichlet boundary conditions
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides);
int_dofs = setdiff (1:space.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - Ai(int_dofs, drchlt_dofs)*u_drchlt;

bpx = struct ('Ai', Ai, 'rhs', rhs, 'int_dofs', int_dofs, 'ndof', space.ndof, 'new_dofs', int_dofs, 'Pi', [], 'Qi', []);
sp{1} = space;

for ilev = 2:nlevels
  msh = msh_refine (msh, 2*ones(1,msh.ndim)); % Only valid in the uniform case
  [space,Proj] = sp_refine (space, msh, 2*ones(1,msh.ndim), degree, regularity); % Uniform refinement
  sp{ilev} = space;

% Assemble the matrices
  Ai = op_gradu_gradv_tp (space, space, msh);
  rhs = op_f_v_tp (space, msh, f);

% Apply Dirichlet boundary conditions
  [u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides);
  int_dofs = setdiff (1:space.ndof, drchlt_dofs);
  rhs(int_dofs) = rhs(int_dofs) - Ai(int_dofs, drchlt_dofs)*u_drchlt;
  u = zeros (space.ndof, 1);
  u(drchlt_dofs) = u_drchlt;

  Caux = 1;
  for idim = 1:msh.ndim
    Caux = kron (Proj{idim}, Caux);
  end
  bpx(ilev-1).Pi = Caux;

  bpx(ilev).Ai = Ai;
  bpx(ilev).rhs = rhs;
  bpx(ilev).int_dofs = int_dofs;
  bpx(ilev).ndof = space.ndof;


  bpx(ilev).Qi = speye (bpx(ilev).ndof);
  for lind = ilev-1:-1:1
    bpx(lind).Qi = bpx(lind+1).Qi * bpx(lind).Pi;
  end

  bpx(ilev).new_dofs = int_dofs; % To be corrected for local refinement

  for lind = 1:ilev
    bpx(lind).Qi = bpx(lind).Qi(:,bpx(lind).new_dofs);
  end


% % Solve the linear system
  int_dofs = bpx(ilev).int_dofs;
  A = bpx(ilev).Ai(int_dofs, int_dofs);
  b = bpx(ilev).rhs(int_dofs);
%   condA(ilev) = eigs(A, 1, 'LM') / eigs(A, 1, 'SM')
  
  tol = 1e-5 / 2^(ilev*degree(1));

  [u(int_dofs), flag, relres, iter, resvec, eigest_j] = ...
    pcg_w_eigest (A, b, tol, numel (int_dofs), @prec_bpx_new, [], bpx, ilev, 1);
%     eigest_j = my_bpx (A, bpx, ilev, 1);
  lambda_min_jac(ilev) = eigest_j(1);
  lambda_max_jac(ilev) = eigest_j(2);
  CondNum_PrecA_jac(ilev) = eigest_j(2) / eigest_j(1)

   [u(int_dofs), flag, relres, iter, resvec, eigest_gs] = ...
     pcg_w_eigest (A, b, tol, numel (int_dofs), @prec_bpx_new, [], bpx, ilev, 2);
 %     eigest_gs = my_bpx (A, bpx, ilev, 2);
   lambda_min_gs(ilev) = eigest_gs(1);
   lambda_max_gs(ilev) = eigest_gs(2);
   CondNum_PrecA_gs(ilev) = eigest_gs(2) / eigest_gs(1)

% %   [u(int_dofs), flag, relres, iter, resvec, eigest_ric] = ...
% %     pcg_w_eigest (A, b, tol, numel (int_dofs), @prec_bpx_new, [], bpx, ilev, 4, msh.ndim);
% %   lambda_min_ric(ilev) = eigest_ric(1);
% %   lambda_max_ric(ilev) = eigest_ric(2);
% %   CondNum_PrecA_ric(ilev) = eigest_ric(2) / eigest_ric(1)

  nint_dofs(ilev) = numel (bpx(ilev).int_dofs);

% %   [error_h1(ilev), error_l2(ilev)] = sp_h1_error (space, msh, u, uex, graduex)

%   save results_bsplines condA CondNum_PrecA_jac CondNum_PrecA_gs nint_dofs lambda_min_jac lambda_max_jac lambda_min_gs lambda_max_gs
% save (filename, 'condA', 'nint_dofs', 'CondNum_PrecA_jac', 'eigest_j', 'iter', 'CondNum_PrecA_jac2', 'eigest_j2', 'iter2')
save (filename, 'nint_dofs', 'CondNum_PrecA_jac', 'lambda_min_jac', 'lambda_max_jac', 'iter')
end

msh_1 = msh; space_1 = space;

end
