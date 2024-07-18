% COARSENING_ALGORITHM: perform a coarsening algorithm based on some indicator,
%  and project the variable in the hierarchical space of the coarsened mesh.
%  Used for Cahn-Hilliard equation, it projects two variables. It also
%  includes a penalization term on the boundary for the projection.
%
% INPUT:
%
%  est:        the error indicator
%  hmsh:       mesh object (see hierarchical_mesh or hierarchical_mesh_mp).
%  hspace:     space object (see hierarchical_space or hierarchical_space_mp_C1).
%  u_n:        field at the previous time step.
%  u_dotn:     time derivative at the previous time step.
%  Cpen:       penalization parameter for the penalty term.
%  old_space:  space from previous iterations. If not changed, some matrices are not recomputed.
%  nmnn_sides: sides where to impose the penalty term.
%
% OUTPUT:
%
%  u_n:       field on the coarsened mesh.
%  u_dotn:    time derivative on the coarsened mesh.
%  hspace:    coarsened space object (see hierarchical_space or hierarchical_space_mp_C1).
%  hmsh:      coarsened mesh object (see hierarchical_mesh or hierarchical_mesh_mp).
%  old_space: space of the previous iteration, updated.
%
% Copyright (C) 2023, 2024 Michele Torre, Rafael Vazquez
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

function [u_n1, udot_n1, hspace, hmsh, old_space] = ...
  coarsening_algorithm(est, hmsh, hspace, adaptivity_data, u_n1, udot_n1, pen_proje, old_space, nmnn_sides)

  %------------------------------------------------------------------  
  % mark
% %%%%%%%%%%%%%%%%%%%%%%%%%% TODO: write the help of the function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  [marked, num_marked_coa] = adaptivity_mark_coarsening_cahn_hilliard (est, hmsh, hspace, adaptivity_data);

  if (num_marked_coa == 0)
    old_space.modified = false;
  else
  %------------------------------------------------------------------
  % coarsening
    hmsh_fine = hmsh;
    hspace_fine = hspace;
    [hmsh, hspace] = adaptivity_coarsen (hmsh, hspace, marked, adaptivity_data);    
        
    if (hspace.ndof == hspace_fine.ndof)
      old_space.modified = false;
    else
      % recompute control variables:  (mass+penalty) \ (G * u)
      [u_n1, udot_n1] = compute_control_variables_coarse_mesh...
        (hmsh, hspace, hmsh_fine, hspace_fine, u_n1, udot_n1, pen_proje, nmnn_sides);

      old_space = struct ('modified', true, 'space', [], 'mesh', [], 'mass_mat', [], ...
                          'lapl_mat', [], 'bnd_mat', [], 'Pen', [], 'pen_rhs', []);
    end
    clear hmsh_fine hspace_fine

  end
end

%--------------------------------------------------------------------------
% compute control variables on the coarser mesh by means of L2-projection
%--------------------------------------------------------------------------
function [u_n_coa, udot_n_coa] = ...
  compute_control_variables_coarse_mesh(hmsh, hspace, hmsh_fine, hspace_fine, u_n, udot_n, pen_proje, nmnn_sides)

  mass_coarse = op_u_v_hier(hspace,hspace,hmsh);

  % penalty term (matrix and vector). The right-hand side is set to zero.
  [Pen, ~] = op_penalty_dudn (hspace, hmsh, nmnn_sides, pen_proje);
  mass_coarse = mass_coarse + Pen;

  hspace_in_hmsh_fine = hspace_in_finer_mesh(hspace, hmsh, hmsh_fine);
  rhs_u = op_Gu_hier (hspace_in_hmsh_fine, hmsh_fine, hspace_fine, u_n);
  u_n_coa = mass_coarse\rhs_u;

  rhs_udot = op_Gu_hier (hspace_in_hmsh_fine, hmsh_fine, hspace_fine, udot_n);
  udot_n_coa = mass_coarse\rhs_udot;

end

%--------------------------------------------------------------------------
% Compute a rhs-vector of the form b_i = (B_i, u) = (B_i, sum_j u_j D_j), 
% with B_i a basis of one (coarse) space, and u written as a combination of 
% basis functions of a second (fine) space, using the finer mesh
%--------------------------------------------------------------------------
function rhs = op_Gu_hier (hspace, hmsh_fine, hspace_fine, uhat_fine)

  rhs = zeros (hspace.ndof, 1);

  ndofs = 0;
  for ilev = 1:hmsh_fine.nlevels
    ndofs = ndofs + hspace.ndof_per_level(ilev);
    if (hmsh_fine.nel_per_level(ilev) > 0)
      x = cell (hmsh_fine.rdim, 1);
      for idim = 1:hmsh_fine.rdim
        x{idim} = reshape (hmsh_fine.msh_lev{ilev}.geo_map(idim,:,:), hmsh_fine.msh_lev{ilev}.nqn, hmsh_fine.nel_per_level(ilev));
      end

      sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh_fine.msh_lev{ilev}, 'value', true);      
      sp_lev = change_connectivity_localized_Csub (sp_lev, hspace, ilev);

      % coefficients
      ndof_until_lev = sum (hspace_fine.ndof_per_level(1:ilev));
      uhat_lev = hspace_fine.Csub{ilev} * uhat_fine(1:ndof_until_lev);
      sp_lev_fine = sp_evaluate_element_list (hspace_fine.space_of_level(ilev), hmsh_fine.msh_lev{ilev}, 'value', true); 
      sp_lev_fine = change_connectivity_localized_Csub (sp_lev_fine, hspace_fine, ilev);
      utemp = sp_eval_msh (uhat_lev, sp_lev_fine, hmsh_fine.msh_lev{ilev}, {'value'});
      u_fine = utemp{1};

      b_lev = op_f_v (sp_lev, hmsh_fine.msh_lev{ilev}, u_fine);

      dofs = 1:ndofs;
      rhs(dofs) = rhs(dofs) + hspace.Csub{ilev}.' * b_lev;
    end
  end

end
