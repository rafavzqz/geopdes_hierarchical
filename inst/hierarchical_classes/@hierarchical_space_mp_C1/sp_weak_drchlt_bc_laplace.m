% SP_WEAK_DRCHLT_BC_LAPLACE: compute the matrix and right hand-side to impose
% the Dirichlet boundary conditions in weak form for Laplace (Poisson) problem.
%
% The code computes the following terms in the left hand-side
% 
%  - \int_{Gamma_D} mu*{du/dn v - dv/dn u + (Cpen /  he) * (u v)}
%
% and in the right hand-side
%
%  - \int_{Gamma_D} mu*{dv/dn g + (Cpen / he) * (v g)}
%
% with u the trial function, v the test function, he the normal characteristic length, 
%  and g the boundary condition to be imposed.
%
%
%   [N_mat, N_rhs] = sp_weak_drchlt_bc_laplace (space, msh, bnd_sides, bnd_func, coeff, Cpen)
%
% INPUTS:
%     
%    space:     space object (see sp_vector)
%    msh:       mesh object (see msh_cartesian)
%    bnd_sides: boundary sides on which the Dirichlet condition is imposed
%    bnd_func:  the condition to be imposed (g in the equations)
%    coeff:     function handle for the viscosity coefficient (mu in the equation)
%    Cpen:      a penalization term
%   
% OUTPUT:
%
%     N_mat:       the computed matrix, to be added in the left hand-side
%     N_rhs:       the computed right hand-side
%
% Copyright (C) 2014 Adriano Cortes, Rafael Vazquez
% Copyright (C) 2015, 2017, 2022 Rafael Vazquez
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

function [A, rhs] = sp_weak_drchlt_bc_laplace (hspace, hmsh, bnd_sides, bnd_func, coeff, Cpen)

  if (nargin < 6 || isempty (Cpen))
    Cpen = 10 * (min (hspace.space_of_level(1).sp_patch{1}.degree) + 1); 
  end

  A = spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);
  rhs = zeros (hspace.ndof, 1);

  boundaries = hmsh.mesh_of_level(1).boundaries;
  Nbnd = cumsum ([0, boundaries.nsides]);
  last_dof = cumsum (hspace.ndof_per_level);
    
  for iref = bnd_sides
    iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    for ilev = 1:hmsh.boundary.nlevels
      patch_bnd_shifting = cumsum ([0 hmsh.boundary.mesh_of_level(ilev).nel_per_patch]);
      patch_shifting = cumsum ([0 hmsh.mesh_of_level(ilev).nel_per_patch]);
      dofs_on_lev = 1:last_dof(ilev);

      for ii = 1:numel(iref_patch_list)
        iptc_bnd = iref_patch_list(ii);
        iptc = boundaries(iref).patches(ii);
        iside = boundaries(iref).faces(ii);
        elems_patch = patch_bnd_shifting(iptc_bnd)+1:patch_bnd_shifting(iptc_bnd+1);
        [~, ~, elements] = intersect (hmsh.boundary.active{ilev}, elems_patch);

        if (~isempty (elements))
          msh_patch = hmsh.mesh_of_level(ilev).msh_patch{iptc};

          msh_side = msh_eval_boundary_side (msh_patch, iside, elements);
          msh_side_from_interior = msh_boundary_side_from_interior (msh_patch, iside);

          sp_bnd = hspace.space_of_level(ilev).sp_patch{iptc}.constructor (msh_side_from_interior);
          msh_side_from_interior_struct = msh_evaluate_element_list (msh_side_from_interior, elements);
          sp_bnd = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', true, 'gradient', true);

          for idim = 1:hmsh.rdim
            x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
          end
          coeff_at_qnodes = coeff (x{:});

          [Cpatch, Cpatch_cols_lev] = sp_compute_Cpatch (hspace.space_of_level(ilev), iptc);
          [~,Csub_rows,Cpatch_cols] = intersect (hspace.Csub_row_indices{ilev}, Cpatch_cols_lev);
          Caux = Cpatch(:,Cpatch_cols) * hspace.Csub{ilev}(Csub_rows,:);
          B = op_gradv_n_u (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
          B = Caux.' * B * Caux;

          g_times_coeff = bnd_func(x{:},iref) .* coeff_at_qnodes;
          gradv_n_g = op_gradv_n_f (sp_bnd, msh_side, g_times_coeff);

% For simplicity I replaced charlen by 2^(-level)
%           coeff_at_qnodes =  coeff_at_qnodes * Cpen ./ msh_side.charlen;
          coeff_at_qnodes =  coeff_at_qnodes * Cpen / 2^(-ilev);
          C = op_u_v (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
          C = Caux.' * C * Caux;

          g_times_coeff = bnd_func(x{:}, iref) .* coeff_at_qnodes;
          g_cdot_v = op_f_v (sp_bnd, msh_side, g_times_coeff);

          A(dofs_on_lev,dofs_on_lev) = A(dofs_on_lev, dofs_on_lev) + (B + B.' - C);
          rhs(dofs_on_lev) = rhs(dofs_on_lev) + Caux.' * (g_cdot_v - gradv_n_g);
        end
      end
    end
  end

% % Compute the matrices to impose the tangential boundary condition weakly
%   for iside = bnd_sides
% 
%     msh_side = msh_eval_boundary_side (msh, iside);
%     msh_side_from_interior = msh_boundary_side_from_interior (msh, iside);
% 
%     sp_bnd = space.constructor (msh_side_from_interior);
%     sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true, 'gradient', true);
% 
%     for idim = 1:msh.rdim
%       x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
%     end
% 
%     coeff_at_qnodes = coeff (x{:});
% 
%     % Since trial and test spaces are the same, we can use B'
%     B = op_gradv_n_u (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
% 
%     g_times_coeff = bnd_func(x{:},iside) .* coeff_at_qnodes;
%     gradv_n_g = op_gradv_n_f (sp_bnd, msh_side, g_times_coeff);
% 
%     coeff_at_qnodes =  coeff_at_qnodes * Cpen ./ msh_side.charlen;
%     C = op_u_v (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
% 
%     g_times_coeff = bnd_func(x{:}, iside) .* coeff_at_qnodes;
%     g_cdot_v = op_f_v (sp_bnd, msh_side, g_times_coeff);
% 
%     A = A + (B + B' - C);
%     rhs = rhs - gradv_n_g + g_cdot_v;
%   end

end
