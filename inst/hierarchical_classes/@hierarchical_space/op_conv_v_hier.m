function rhs = op_conv_v_hier (hspace, hmsh, conv, u)

  rhs = zeros (hspace.ndof, 1);

  ndofs = 0;
  for ilev = 1:hmsh.nlevels
    ndofs = ndofs + hspace.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)
      x = cell (hmsh.rdim, 1);
      for idim = 1:hmsh.rdim
        x{idim} = reshape (hmsh.msh_lev{ilev}.geo_map(idim,:,:), hmsh.mesh_of_level(ilev).nqn, hmsh.nel_per_level(ilev));
      end
      sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
      [eu, ~] = hspace_eval_hmsh(u, hspace, hmsh);
      
      b_lev = op_f_v (sp_lev, hmsh.msh_lev{ilev}, conv(x{:},eu(:,1:hmsh.msh_lev{ilev}.nel)));

      dofs = 1:ndofs;
      rhs(dofs) = rhs(dofs) + hspace.Csub{ilev}.' * b_lev;
    end
  end

end

