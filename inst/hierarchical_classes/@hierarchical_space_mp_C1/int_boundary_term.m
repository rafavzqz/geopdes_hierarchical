%--------------------------------------------------------------------------
% Boundary term, \int_\Gamma (\Delta u) (\partial v / \partial n)
%--------------------------------------------------------------------------
function [A] = int_boundary_term (hspace, hmsh, lambda, nmnn_sides)

  if (~isempty(nmnn_sides))

    A = spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);

    boundaries = hmsh.mesh_of_level(1).boundaries;
    Nbnd = cumsum ([0, boundaries.nsides]);
    last_dof = cumsum (hspace.ndof_per_level);

    for iref = nmnn_sides
      iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
      for ilev = 1:hmsh.boundary.nlevels
        patch_bnd_shifting = cumsum ([0 hmsh.boundary.mesh_of_level(ilev).nel_per_patch]);
%        patch_shifting = cumsum ([0 hmsh.mesh_of_level(ilev).nel_per_patch]);
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
            sp_bnd = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', false, 'gradient', true, 'laplacian', true);

            for idim = 1:hmsh.rdim
              x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
            end
            coe_side = lambda (x{:});

            [Cpatch, Cpatch_cols_lev] = sp_compute_Cpatch (hspace.space_of_level(ilev), iptc);
            [~,Csub_rows,Cpatch_cols] = intersect (hspace.Csub_row_indices{ilev}, Cpatch_cols_lev);
            Caux = Cpatch(:,Cpatch_cols) * hspace.Csub{ilev}(Csub_rows,:);

            tmp = op_gradv_n_laplaceu(sp_bnd ,sp_bnd ,msh_side, coe_side);

            A(dofs_on_lev,dofs_on_lev) = A(dofs_on_lev, dofs_on_lev) + Caux.' * tmp * Caux;
          end
        end
      end
    end

  else
    A = [];
  end
end
