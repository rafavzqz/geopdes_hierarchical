%--------------------------------------------------------------------------
% penalty term
%--------------------------------------------------------------------------
function [P, rhs] = penalty_matrix (hspace, hmsh, nmnn_sides, pen)
    
  P =  spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);
  rhs = zeros(hspace.ndof,1);

  for iside = 1:length(nmnn_sides)
    [mass_pen, rhs_pen] = penalty_grad (hspace, hmsh, nmnn_sides(iside), pen);
    P = P + mass_pen; 
    rhs = rhs + rhs_pen;
  end

end

%--------------------------------------------------------------------------
function [mass_pen, rhs_pen] = penalty_grad (hspace, hmsh, nmnn_side, pen)

  hmsh_side_int = hmsh_boundary_side_from_interior (hmsh, nmnn_side ) ;
  shifting_indices = cumsum ([0 hmsh.boundary(nmnn_side).nel_per_level]);

  mass_pen = spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);

  ndofs_u = 0;
  ndofs_v = 0;
  for ilev = 1:hmsh.boundary(nmnn_side).nlevels

    ndofs_u = ndofs_u + hspace.ndof_per_level(ilev);
    ndofs_v = ndofs_v + hspace.ndof_per_level(ilev);
    
    if (hmsh.boundary(nmnn_side).nel_per_level(ilev) > 0)

      % mesh of the selected side
      elements = shifting_indices(ilev)+1:shifting_indices(ilev+1); 
      hmsh_side = hmsh_eval_boundary_side (hmsh, nmnn_side, elements);

      % reconstruct the part of the space defined on the selected boundary
      msh_side_from_interior = hmsh_side_int.mesh_of_level(ilev);
      sp_bnd = hspace.space_of_level(ilev).constructor (msh_side_from_interior);
      msh_side_from_interior_struct = msh_evaluate_element_list (msh_side_from_interior, hmsh_side_int.active{ilev});

      % evaluate the space
      spu_lev = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'gradient', true);
      spu_lev = change_connectivity_localized_Csub (spu_lev, hspace, ilev);

      % compute the matrix
      coe_side = repmat(hmsh_side.element_size, hmsh_side.nqn,1);
      coe_side = pen .* coe_side;
      K_lev = op_gradu_n_gradv_n(spu_lev, spu_lev, hmsh_side, coe_side);

      dofs_u = 1:ndofs_u;
      dofs_v = 1:ndofs_v;
      Ktmp =  hspace.Csub{ilev}.' * K_lev * hspace.Csub{ilev};
      mass_pen(dofs_v,dofs_u) = mass_pen(dofs_v,dofs_u) + Ktmp;
    end
  end

  rhs_pen  = zeros(hspace.ndof,1);

end

