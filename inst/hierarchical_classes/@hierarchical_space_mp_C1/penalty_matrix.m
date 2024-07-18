%--------------------------------------------------------------------------
% penalty term
%--------------------------------------------------------------------------
function [A, rhs] = penalty_matrix (hspace, hmsh, nmnn_sides, Cpen)

  A = spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);
  rhs = zeros (hspace.ndof, 1);

  boundaries = hmsh.mesh_of_level(1).boundaries;
  Nbnd = cumsum ([0, boundaries.nsides]);
  last_dof = cumsum (hspace.ndof_per_level);

  for iref = nmnn_sides
    iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    for ilev = 1:hmsh.boundary.nlevels
      patch_bnd_shifting = cumsum ([0 hmsh.boundary.mesh_of_level(ilev).nel_per_patch]);
%      patch_shifting = cumsum ([0 hmsh.mesh_of_level(ilev).nel_per_patch]);
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

          % For simplicity I replaced charlen by 2^(-level)
%           coeff_at_qnodes =  coeff_at_qnodes * Cpen ./ msh_side.charlen;
          coeff_at_qnodes = Cpen * ones(msh_side.nqn, msh_side.nel)/ 2^(-ilev);

          [Cpatch, Cpatch_cols_lev] = sp_compute_Cpatch (hspace.space_of_level(ilev), iptc);
          [~,Csub_rows,Cpatch_cols] = intersect (hspace.Csub_row_indices{ilev}, Cpatch_cols_lev);
          Caux = Cpatch(:,Cpatch_cols) * hspace.Csub{ilev}(Csub_rows,:);

          tmp = op_gradu_n_gradv_n(sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
          tmp = Caux.' * tmp * Caux;

          A(dofs_on_lev,dofs_on_lev) = A(dofs_on_lev, dofs_on_lev) + tmp;
        end
      end
    end
  end

end

