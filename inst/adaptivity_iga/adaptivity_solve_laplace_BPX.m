% ADAPTIVITY_SOLVE_LAPLACE: assemble and solve the linear system for Laplacian problem, using hierarchical spaces.
%
% The function solves the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^n)
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% USAGE:
%
% u = adaptivity_solve_laplace (hmsh, hspace, method_data)
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

function [u, bpx] = adaptivity_solve_laplace_gerarchia (hmsh, hspace, problem_data, bpx)

stiff_mat = op_gradu_gradv_hier (hspace, hspace, hmsh, problem_data.c_diff);
rhs = op_f_v_hier (hspace, hmsh, problem_data.f);
% mass = op_u_v_hier (hspace, hspace, hmsh, problem_data.c_diff);

% Apply Neumann boundary conditions
if (~isfield (struct (hmsh), 'npatch')) % Single patch case
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

% Apply Dirichlet boundary conditions
u = zeros (hspace.ndof, 1);
[u_dirichlet, dirichlet_dofs] = sp_drchlt_l2_proj (hspace, hmsh, problem_data.h, problem_data.drchlt_sides);
u(dirichlet_dofs) = u_dirichlet;

int_dofs = setdiff (1:hspace.ndof, dirichlet_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, dirichlet_dofs)*u(dirichlet_dofs);

% Initialize
  hmsh_bpx(1)     = hierarchical_mesh (hmsh.mesh_of_level(0), method_data.nsub_refine); % Bisogna passare method_data
  hspace_bpx(1)   = hierarchical_space (hmsh_bps, hspace.space_of_level(0), method_data.space_type, method_data.truncated);

% CI MANCA L'INFORMAZIONE DI METHOD_DATA.BPX_DOFS  
for ref = 1:hmsh.nlevels-1
  marked = cell(ref,1);
  marked{ref} = hmsh.deactivated{ref};
  adap_data.flag = 'elements';
  [hmsh_bpx(ref+1), hspace_bpx(ref+1), Cref, new_dofsXXX] = adaptivity_refine (hmsh, hspace, marked, adap_flag); % raffinando per elementi
  
  bpx(ref).Pi = Cref;
%% Calcolare queste bene  
%  bpx(ref+1).Ai = stiff_mat;
%  bpx(ref+1).rhs = rhs;
%  bpx(ref+1).int_dofs = int_dofs;
  
% Aggiungere calcolo della Qi
  bpx(ref+1).ndof = hspace.ndof;
  if (strcmpi (method_data.bpx_dofs, 'All_dofs'))
    bpx(ref+1).new_dofs = int_dofs;
  else
    bpx(ref+1).new_dofs = intersect (int_dofs, hspace.active{ref+1}XXXXXX); % Deve essere numerazione di hspace
  end

% Controllare che sia fatto nel modo giusto (togliere condizioni di bordo)  
  for lind = 1:ref+1
    bpx(lind).Qi = bpx(lind).Qi(:,bpx(lind).new_dofs);
  end

end



% Solve the linear system
% u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);
  tol = 1e-10;% / 2^(ilev*degree(1));
  
  A = bpx(end).Ai(int_dofs, int_dofs);
  b = bpx(end).rhs(int_dofs);

  [u(int_dofs), flag, relres, iter, resvec, eigest_j] = ...
    pcg_w_eigest (A, b, tol, numel (int_dofs), @prec_bpx_new, [], bpx, this_level, 1);
%     eigest_j = my_bpx (A, bpx, ilev, 1);
  lambda_min_jac = eigest_j(1);
  lambda_max_jac = eigest_j(2);
  CondNum_PrecA_jac = eigest_j(2) / eigest_j(1);
