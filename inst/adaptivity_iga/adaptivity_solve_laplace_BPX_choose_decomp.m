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

function [u, bpx, solution_data] = adaptivity_solve_laplace_BPX_choose_decomp (hmsh, hspace, problem_data, method_data, decomp, solution_data)

geometry = geo_load (problem_data.geo_name);

stiff_mat = op_gradu_gradv_hier_lowmem (hspace, hspace, hmsh, problem_data.c_diff);
rhs = op_f_v_hier (hspace, hmsh, problem_data.f);

mass_mat = op_u_v_hier_lowmem (hspace, hspace, hmsh);

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
msh_one = hmsh.mesh_of_level(1);
msh_one = msh_cartesian (msh_one.breaks, msh_one.qn, msh_one.qw, geometry);

% hmsh_bpx   = hierarchical_mesh (hmsh.mesh_of_level(1), method_data.nsub_refine); % Bisogna passare method_data
hmsh_bpx   = hierarchical_mesh (msh_one, method_data.nsub_refine); % Bisogna passare method_data
hspace_bpx = hierarchical_space (hmsh_bpx, hspace.space_of_level(1), method_data.space_type, method_data.truncated);
bpx(1).Qi = speye (hspace_bpx(1).ndof);
Cref = speye (hspace_bpx(1).ndof);

for ref = 1:hmsh.nlevels-1
  bpx(ref).Ai = op_gradu_gradv_hier_lowmem (hspace_bpx(ref), hspace_bpx(ref), hmsh_bpx(ref));
  bpx(ref).rhs = op_f_v_hier (hspace_bpx(ref), hmsh_bpx(ref), problem_data.f);
  bpx(ref).ndof = hspace_bpx(ref).ndof;
  bpx(ref).Mi = op_u_v_hier_lowmem (hspace_bpx(ref), hspace_bpx(ref), hmsh_bpx(ref));
  
  cumndof = cumsum ([0 hspace_bpx(ref).ndof_per_level]);
  [~, drchlt_dofs_lev] = sp_drchlt_l2_proj (hspace_bpx(ref), hmsh_bpx(ref), problem_data.h, problem_data.drchlt_sides);
  bpx(ref).int_dofs = setdiff (1:hspace_bpx(ref).ndof, drchlt_dofs_lev);

% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Here we decide if we take all the functions until level ref (All_dofs)
%  only the newly added functions until level ref (New_dofs)
%  or the functions newly added + truncated (Mod_dofs)
% In the last two cases, we have to change the variable new_dofs that
%  is given by adaptivity_refine
% Forse bpx_dofs va messo su adaptivity_data (e adap_data), e non su method_data.
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  for ide = 1:numel(decomp)
    if (strcmpi (decomp{ide}, 'All_dofs'))
      bpx(ref).new_dofs{ide} = bpx(ref).int_dofs;
    elseif (strcmpi (decomp{ide}, 'New_dofs')) % XXXXXXX This is only valid for the non-truncated basis
      new_dofs = cumndof(ref) + (1:numel(hspace_bpx(ref).active{end}));
      bpx(ref).new_dofs{ide} = intersect (bpx(ref).int_dofs, new_dofs);
    elseif (strcmpi (decomp{ide}, 'Mod_dofs'))
      new_dofs = cumndof(ref) + (1:numel(hspace_bpx(ref).active{end}));
      if (ref > 1 && hspace.truncated)
        cumndof_old = cumsum ([0 hspace_bpx(ref-1).ndof_per_level]);
        last_level_dofs = (cumndof(ref)+1):cumndof(ref+1);
      
        truncated_dofs = [];
        for lev = 1:hmsh_bpx(ref-1).nlevels
          [~,ind_new,ind_old] = intersect (hspace_bpx(ref).active{lev}, hspace_bpx(ref-1).active{lev});
          global_indices_old = cumndof_old(lev) + ind_old;
          global_indices = cumndof(lev) + ind_new;
          for jj = 1:numel(ind_old)
            [aa,bb] = find(Cref(:,global_indices_old(jj)));
            bb = intersect (aa, last_level_dofs);
            if (~isempty (bb))
              truncated_dofs = union (truncated_dofs, global_indices(jj));
            end
          end
        end
      
        new_dofs = union (new_dofs, truncated_dofs);
      end
%     new_dofs = union (new_dofs, new_dofs-method_data.degree(1));
      bpx(ref).new_dofs{ide} = intersect (bpx(ref).int_dofs, new_dofs);
    
    elseif (strcmpi (decomp{ide}, 'Support_dofs'))
% Functions of previous levels entering in \Omega_{l+1}
      new_cells = hmsh_bpx(ref).active{end};
      funs = sp_get_basis_functions (hspace_bpx(ref).space_of_level(end), hmsh_bpx(ref).mesh_of_level(end), new_cells);
      [~,jj] = find (hspace_bpx(ref).Csub{ref}(funs,:)); 
      new_dofs = unique(jj);
      bpx(ref).new_dofs{ide} = intersect (bpx(ref).int_dofs, new_dofs);
    end
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

% Check that it is done correctly (remove boundary conditions)
bpx(hmsh.nlevels).Ai = stiff_mat;
bpx(hmsh.nlevels).rhs = rhs;
bpx(hmsh.nlevels).Pi = [];
bpx(hmsh.nlevels).ndof = hspace.ndof;
bpx(hmsh.nlevels).int_dofs = int_dofs;
bpx(hmsh.nlevels).Mi = mass_mat;

for ide = 1:numel(decomp)
  if (strcmpi (decomp{ide}, 'All_dofs'))
    bpx(hmsh.nlevels).new_dofs{ide} = int_dofs;
  elseif (strcmpi (decomp{ide}, 'New_dofs'))
    cumndof = cumsum ([0 hspace.ndof_per_level]);
    new_dofs = cumndof(hmsh.nlevels) + (1:numel(hspace.active{end}));
    bpx(hmsh.nlevels).new_dofs{ide} = intersect (int_dofs, new_dofs);
  elseif (strcmpi (decomp{ide}, 'Mod_dofs')) % XXXXXXXX Only for this particular 1D case
    cumndof = cumsum ([0 hspace.ndof_per_level]);
    new_dofs = cumndof(hmsh.nlevels) + (1:numel(hspace.active{end}));
    if (hspace.nlevels > 1 && hspace.truncated)
      cumndof_old = cumsum ([0 hspace_bpx(end-1).ndof_per_level]);
      last_level_dofs = (cumndof(end-1)+1):cumndof(end);
      
      truncated_dofs = [];
      for lev = 1:hmsh_bpx(end-1).nlevels
        [~,ind_new,ind_old] = intersect (hspace_bpx(end).active{lev}, hspace_bpx(end-1).active{lev});
        global_indices_old = cumndof_old(lev) + ind_old;
        global_indices = cumndof(lev) + ind_new;
        for jj = 1:numel(ind_old)
          [aa,bb] = find(Cref(:,global_indices_old(jj)));
          bb = intersect (aa, last_level_dofs);
          if (~isempty (bb))
            truncated_dofs = union (truncated_dofs, global_indices(jj));
          end
        end
      end
      
      new_dofs = union (new_dofs, truncated_dofs);
    end
    bpx(hmsh.nlevels).new_dofs{ide} = intersect (int_dofs, new_dofs);

  elseif (strcmpi (decomp{ide}, 'Support_dofs'))
% Functions of previous levels entering in \Omega_{l+1}
    new_cells = hmsh.active{end};
    funs = sp_get_basis_functions (hspace.space_of_level(end), hmsh.mesh_of_level(end), new_cells);
    [~,jj] = find (hspace.Csub{end}(funs,:)); 
    new_dofs = unique(jj);
    bpx(hmsh.nlevels).new_dofs{ide} = intersect (bpx(hmsh.nlevels).int_dofs, new_dofs);
  end
end
bpx(hmsh.nlevels).Qi = speye (hspace.ndof);

for lind = 1:hmsh.nlevels
  Qi{lind} = bpx(lind).Qi;
  new_dofs_cell{lind} = bpx(lind).new_dofs;
end


% Solve the linear system
% u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);

if (~isfield (method_data, 'cond_iter') || strcmpi (method_data.cond_iter, 'cond'))
  tol = method_data.tol / 2^(hmsh.nlevels);%*method_data.degree(1));
elseif (strcmpi (method_data.cond_iter, 'iter'))
  tol = method_data.tol;
end
  
A = bpx(end).Ai(int_dofs, int_dofs);
b = bpx(end).rhs(int_dofs);


for ide = 1:numel(decomp)
  
  for lind = 1:hmsh.nlevels
    bpx(lind).new_dofs = new_dofs_cell{lind}{ide};
    bpx(lind).Qi = Qi{lind}(:,bpx(lind).new_dofs);
  end
  
% % 1 is Jacobi, 2 is Gauss-Seidel, 4 is Richardson
% % JACOBI SMOOTHER
%   [u(int_dofs), flag, relres, iter, resvec, eigest_j] = ...
%     pcg_w_eigest (A, b, tol, numel (int_dofs), @prec_bpx_new, [], bpx, hmsh.nlevels, 1);
%   CondNum_PrecA_jac = eigest_j(2) / eigest_j(1);
%   disp (['Condition number for BPX with Jacobi: ', num2str(CondNum_PrecA_jac)]);
%   disp(['Eigenvalues: ', num2str(eigest_j)])
%   disp(['Iterations: ', num2str(iter)])
%   niter_jac = iter;
% %   CondNum_PrecA_jac = 0;

% Remove fine empty levels (elements without active functions)
  finest_level = find(hspace.ndof_per_level>0,1,'last');
  bpx = bpx(1:finest_level);
  
  eigest_j = [0 0]; CondNum_PrecA_jac = 0; niter_jac = 0;

% % 1 is Jacobi, 2 is Gauss-Seidel, 4 is Richardson
% % GAUSS-SEIDEL SMOOTHER
tt = cputime;
  [u(int_dofs), flag, relres, iter, resvec, eigest_gs] = ...
    pcg_w_eigest (A, b, tol, numel (int_dofs), @prec_bpx_new, [], bpx, finest_level, 2);
ttime = cputime - tt;
  CondNum_PrecA_gs = eigest_gs(2) / eigest_gs(1);
  disp (['Condition number for BPX with Gauss-Seidel: ', num2str(CondNum_PrecA_gs)]);
  disp(['Eigenvalues: ', num2str(eigest_gs)])
  disp(['Iterations: ', num2str(iter)])
  niter_gs = iter;
  %   CondNum_PrecA_gs = 0;
  others = [];% others = 0;
  
  ndof_finest = numel (bpx(end).new_dofs);
  ndof_biggest = max (arrayfun (@(x) numel(x.new_dofs), bpx));
  ndof_average = sum (arrayfun (@(x) numel(x.new_dofs), bpx)) / numel(bpx);
  ndof_levels = arrayfun (@(x) numel(x.new_dofs), bpx);

  solution_data.dec(ide).name = decomp{ide};
  if (hmsh.nlevels == 1)
%     solution_data.dec(ide).Cond_BPX_jac(iter) = CondNum_PrecA_jac;
    solution_data.dec(ide).Cond_BPX_gs = CondNum_PrecA_gs;
%     solution_data.dec(ide).eig_jac = eigest_jac;
    solution_data.dec(ide).eig_gs = eigest_gs;
%     solution_data.dec(ide).niter_jac = niter_jac;
    solution_data.dec(ide).niter_gs = niter_gs;
    solution_data.dec(ide).ndof_finest = ndof_finest;
    solution_data.dec(ide).ndof_biggest = ndof_biggest;
    solution_data.dec(ide).ndof_average = ndof_average;
    solution_data.dec(ide).ndof_levels{1} = ndof_levels;
    solution_data.dec(ide).cpu_time = ttime;
  else
%     solution_data.dec(ide).Cond_BPX_jac(end+1) = CondNum_PrecA_jac;
    solution_data.dec(ide).Cond_BPX_gs(end+1) = CondNum_PrecA_gs;
%     solution_data.dec(ide).eig_jac(end+1,:) = eigest_jac;
    solution_data.dec(ide).eig_gs(end+1,:) = eigest_gs;
%     solution_data.dec(ide).niter_jac(end+1,:) = niter_jac;
    solution_data.dec(ide).niter_gs(end+1,:) = niter_gs;
    solution_data.dec(ide).ndof_finest(end+1,:) = ndof_finest;
    solution_data.dec(ide).ndof_biggest(end+1,:) = ndof_biggest;
    solution_data.dec(ide).ndof_average(end+1,:) = ndof_average;
    solution_data.dec(ide).ndof_levels{end+1} = ndof_levels;
    solution_data.dec(ide).cpu_time(end+1,:) = ttime;
  end
  
end

% CondA = 0;
%   [~,eigest_A(1),flag1] = eigs(A,1,'SM');
%   [~,eigest_A(2),flag2] = eigs(A,1,'LM');
%   eigest_A = [eigs(A,1,'SM'), eigs(A,1,'LM')];
%   if (flag1 || flag2)
  [~, flag, relres, niter_A, resvec, eigest_A] = ...
    pcg_w_eigest (A, b, tol, numel (int_dofs));
    CondA = eigest_A(2) / eigest_A(1)
    
  if (hmsh.nlevels == 1)
    solution_data.dec(ide).CondA = CondA;
    solution_data.dec(ide).eig_A = eigest_A;
    solution_data.dec(ide).niter_A = niter_A;
  else
    solution_data.dec(ide).CondA(end+1) = CondA;
    solution_data.dec(ide).eig_A(end+1,:) = eigest_A;
    solution_data.dec(ide).niter_A(end+1) = niter_A;
  end
  
  u(int_dofs) = A \ b;
%   else
% eigest_A = [0 0];
%     CondA = cond(full(A))
%   end
%   eigv = eig(full(A));
%   CondA = max(eigv)/min(eigv)
