function A = int_boundary_term (hspace, hmsh, lambda, nmnn_sides)
  A =  spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);    
  for iside = 1:length(nmnn_sides)   
    A = A + op_gradv_n_laplaceu_hier(hspace, hmsh, nmnn_sides(iside), lambda );        
  end
end

%--------------------------------------------------------------------------
% OP_GRADVN_LAPLACEU_HIER: assemble the matrix A = [a(i,j)], a(i,j) = (epsilon (grad v n)_j, Delta u_i), 
%  with n the normal vector to the boundary.
%--------------------------------------------------------------------------
function varargout = op_gradv_n_laplaceu_hier (hspace, hmsh, nmnn_side, coeff)

  if (nargin == 3)
    coeff = @(varargin) ones (size(varargin{1}));
  end
    
  % side with info from the interior for space reconstruction
  hmsh_side_int = hmsh_boundary_side_from_interior (hmsh, nmnn_side);
  shifting_indices = cumsum ([0 hmsh.boundary(nmnn_side).nel_per_level]);
    
  K = spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);

  ndofs_u = 0;
  ndofs_v = 0;
  for ilev = 1:hmsh.boundary(nmnn_side).nlevels

    ndofs_u = ndofs_u + hspace.ndof_per_level(ilev);
    ndofs_v = ndofs_v + hspace.ndof_per_level(ilev);

    if (hmsh.boundary(nmnn_side).nel_per_level(ilev) > 0)
      % mesh of the selected side
      elements = shifting_indices(ilev)+1:shifting_indices(ilev+1); 
      hmsh_side = hmsh_eval_boundary_side (hmsh, nmnn_side, elements);

      x = cell (hmsh.rdim,1);
      for idim = 1:hmsh.rdim
        x{idim} = reshape (hmsh_side.geo_map(idim,:,:), hmsh_side.nqn, hmsh_side.nel);
      end
        
      % reconstruct the part of the space defined on the selected boundary
      msh_side_from_interior = hmsh_side_int.mesh_of_level(ilev);
      sp_bnd = hspace.space_of_level(ilev).constructor (msh_side_from_interior);
      msh_side_from_interior_struct = msh_evaluate_element_list (msh_side_from_interior, hmsh_side_int.active{ilev});
        
      % evaluate the space
      spu_lev = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'laplacian', true);
      spv_lev = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', false, 'gradient', true);
        
      spu_lev = change_connectivity_localized_Csub (spu_lev, hspace, ilev);
      spv_lev = change_connectivity_localized_Csub (spv_lev, hspace, ilev);     

      % compute the matrix
      K_lev = op_gradv_n_laplaceu (spu_lev, spv_lev, hmsh_side, coeff (x{:}));

      dofs_u = 1:ndofs_u;
      dofs_v = 1:ndofs_v;
      Ktmp =  hspace.Csub{ilev}.' * K_lev * hspace.Csub{ilev};
      K(dofs_v,dofs_u) = K(dofs_v,dofs_u) +Ktmp;    
    end
  end

  if (nargout == 1)
    varargout{1} = K;
  elseif (nargout == 3)
    [rows, cols, vals] = find (K);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end

