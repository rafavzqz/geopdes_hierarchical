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
% BPX NOTATION
% Ai:       stiffness matrix of level
% rhs:      right-hand side of level
% int_dofs: internal dofs (usually int_dofs)
% ndof:     total number of dofs (usually ndof)
% new_dofs: basis functions (in int_dofs) that were not present in the previous level
% Pi:       projector between two consecutive levels (from l-1 to l)
% Qi:       projector between the current level and the finest one (from l to nlevels)
%
% Copyright (C) 2015, 2016, 2017 Eduardo M. Garau, Rafael Vazquez
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

function [u, bpx] = adaptivity_solve_laplace_BPX2 (hmsh, hspace, problem_data, method_data)

stiff_mat = op_gradu_gradv_hier (hspace, hspace, hmsh, problem_data.c_diff);
rhs = op_f_v_hier (hspace, hmsh, problem_data.f);

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    FROM HERE WE DO THE BPX PRECONDITIONER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
hmsh_bpx   = hierarchical_mesh (hmsh.mesh_of_level(1), method_data.nsub_refine); % Bisogna passare method_data
hspace_bpx = hierarchical_space (hmsh_bpx, hspace.space_of_level(1), method_data.space_type, method_data.truncated);
bpx(1).Qi = speye (hspace_bpx(1).ndof);

% CI MANCA L'INFORMAZIONE DI METHOD_DATA.BPX_DOFS  
for ref = 1:hmsh.nlevels-1
  bpx(ref).Ai = op_gradu_gradv_hier (hspace_bpx(ref), hspace_bpx(ref), hmsh_bpx(ref));
  bpx(ref).rhs = op_f_v_hier (hspace_bpx(ref), hmsh_bpx(ref), problem_data.f);
  bpx(ref).ndof = hspace_bpx(ref).ndof;
  
  drchlt_dofs_lev = [];
  cumndof = cumsum ([0 hspace_bpx(ref).ndof_per_level]);
  for ii = problem_data.drchlt_sides(:)'
%     aux_dofs = hspace_bpx(ref+1).boundary(ii).dofs; % QUESTI NON ESISTONO
    aux_dofs = intersect (hspace.space_of_level(ref).boundary(ii).dofs, hspace_bpx(ref).active{end});
    aux_dofs = aux_dofs + cumndof(ref);
    drchlt_dofs_lev = union (drchlt_dofs_lev, aux_dofs);
  end
  bpx(ref).int_dofs = setdiff (1:hspace_bpx(ref).ndof, drchlt_dofs_lev);

  if (strcmpi (method_data.bpx_dofs, 'All_dofs'))
    bpx(ref).new_dofs = bpx(ref).int_dofs;
  elseif (strcmpi (method_data.bpx_dofs, 'New_dofs')) % XXXXXXX This is only valid for the non-truncated basis
    new_dofs = cumndof(ref) + (1:numel(hspace_bpx(ref).active{end}));
    bpx(ref).new_dofs = intersect (bpx(ref).int_dofs, new_dofs);
  end
  
  marked = cell(ref,1);
  marked{ref} = hmsh.deactivated{ref};
  adap_data.flag = 'elements';
%   [hmsh_bpx(ref+1), hspace_bpx(ref+1), Cref, new_dofsXXX] = adaptivity_refine (hmsh, hspace, marked, adap_data); % raffinando per elementi
  [hmsh_bpx(ref+1), hspace_bpx(ref+1), Cref] = adaptivity_refine (hmsh_bpx(ref), hspace_bpx(ref), marked, adap_data); % raffinando per elementi
  bpx(ref+1).ndof = hspace_bpx(ref+1).ndof;
  
  bpx(ref).Pi = Cref;
  bpx(ref+1).Qi = speye (hspace_bpx(ref+1).ndof);
  for lind = ref:-1:1
    bpx(lind).Qi = bpx(lind+1).Qi * bpx(lind).Pi;
  end
  
end

% Controllare che sia fatto nel modo giusto (togliere condizioni di bordo)  
bpx(hmsh.nlevels).Ai = stiff_mat;
bpx(hmsh.nlevels).rhs = rhs;
bpx(hmsh.nlevels).Pi = [];
bpx(hmsh.nlevels).ndof = hspace.ndof;
bpx(hmsh.nlevels).int_dofs = int_dofs;
if (strcmpi (method_data.bpx_dofs, 'All_dofs'))
  bpx(hmsh.nlevels).new_dofs = int_dofs;
elseif (strcmpi (method_data.bpx_dofs, 'New_dofs'))
  cumndof = cumsum ([0 hspace.ndof_per_level]);
  new_dofs = cumndof(hmsh.nlevels) + (1:numel(hspace.active{end}));
  bpx(hmsh.nlevels).new_dofs = intersect (int_dofs, new_dofs);
end
bpx(hmsh.nlevels).Qi = speye (hspace.ndof);

for lind = 1:hmsh.nlevels
  bpx(lind).Qi = bpx(lind).Qi(:,bpx(lind).new_dofs);
end


% Solve the linear system
% u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);
  tol = 1e-10;% / 2^(ilev*degree(1));
  
  A = bpx(end).Ai(int_dofs, int_dofs);
  b = bpx(end).rhs(int_dofs);

  [u(int_dofs), flag, relres, iter, resvec, eigest_j] = ...
    pcg_w_eigest (A, b, tol, numel (int_dofs), @prec_bpx_new, [], bpx, hmsh.nlevels, 1);
%     eigest_j = my_bpx (A, bpx, ilev, 1);
  lambda_min_jac = eigest_j(1);
  lambda_max_jac = eigest_j(2);
  CondNum_PrecA_jac = eigest_j(2) / eigest_j(1)
