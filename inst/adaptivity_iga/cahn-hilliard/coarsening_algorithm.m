%--------------------------------------------------------------------------
% coarsening
%--------------------------------------------------------------------------
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO: identical to single-patch
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
