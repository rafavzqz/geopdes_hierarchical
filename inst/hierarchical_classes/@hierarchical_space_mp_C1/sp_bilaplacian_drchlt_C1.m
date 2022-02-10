function [u_drchlt, drchlt_dofs] = sp_bilaplacian_drchlt_C1 (hspace, hmsh, refs, h, dudn)

% refs should be the whole boundary, for now
M = spalloc (hspace.ndof, hspace.ndof, hspace.ndof);
rhs = zeros (hspace.ndof, 1);

M2 = spalloc (hspace.ndof, hspace.ndof, hspace.ndof);
rhs2 = zeros (hspace.ndof, 1);

drchlt_dofs = [];
drchlt_dofs2 = [];

boundaries = hmsh.mesh_of_level(1).boundaries;
for iref = refs
  href = @(varargin) h(varargin{:}, iref);
  for bnd_side = 1:boundaries(iref).nsides
    iptc_bnd = sum([boundaries(1:iref-1).nsides]) + bnd_side;
    iptc = boundaries(iref).patches(bnd_side);
    iside = boundaries(iref).faces(bnd_side);

    ndofs = 0;
    for ilev = 1:hmsh.boundary.nlevels
      ndofs = ndofs + hspace.ndof_per_level(ilev);
      active_cells = hmsh.boundary.active{ilev};
      Nind = cumsum ([0 hmsh.boundary.mesh_of_level(ilev).nel_per_patch]);
      [~,active_cells_patch,~] = intersect (Nind(iptc_bnd)+1:Nind(iptc_bnd+1), active_cells);
      
      if (~isempty (active_cells_patch))
        msh_lev_patch = hmsh.mesh_of_level(ilev).msh_patch{iptc};
        sp_lev_patch = hspace.space_of_level(ilev).sp_patch{iptc};
%         msh_bnd = msh_lev_patch.boundary(iside);
        msh_side = msh_eval_boundary_side (msh_lev_patch, iside, active_cells_patch);
        msh_side_from_interior = msh_boundary_side_from_interior (msh_lev_patch, iside);
%       hmsh_sfi = hmsh_boundary_side_from_interior (hmsh, iside);

        sp_bnd = sp_lev_patch.constructor (msh_side_from_interior);
        msh_bnd_struct = msh_evaluate_element_list (msh_side_from_interior, active_cells_patch);
        sp_bnd_struct = sp_evaluate_element_list (sp_bnd, msh_bnd_struct, 'value', true, 'gradient', true);
%         sp_bnd_struct = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true, 'gradient', true);

        sp_bnd = sp_lev_patch.boundary(iside);
%         Cpatch = hspace.space_of_level(ilev).Cpatch{iptc};
        Cpatch = hspace.space_of_level(ilev).Cpatch{iptc}(:,hspace.Csub_row_indices{ilev});
        CC = Cpatch * hspace.Csub{ilev}; % This is probably very unefficient. Who cares...

        [~,icol] = find (CC(sp_bnd.dofs,:));
        [~,jcol] = find (CC(sp_bnd.adjacent_dofs,:));
    
        drchlt_dofs = union (drchlt_dofs, icol);
        drchlt_dofs2 = union (drchlt_dofs2, jcol);
    
        for idim = 1:hmsh.rdim
          x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
        end
        coeff_at_qnodes = ones (size(x{1}));
        dudn_at_qnodes = dudn (x{:},iref);

%        M(1:ndofs,1:ndofs) = M(1:ndofs,1:ndofs) + CC(sp_bnd.dofs,:).' * op_u_v_tp (sp_bnd, sp_bnd, msh_bnd, coeff_at_qnodes) * CC(sp_bnd.dofs,:);
        M(1:ndofs,1:ndofs) = M(1:ndofs,1:ndofs) + CC.' * op_u_v (sp_bnd_struct, sp_bnd_struct, msh_side, coeff_at_qnodes) * CC;
        rhs(1:ndofs) = rhs(1:ndofs) + CC.' * op_f_v (sp_bnd_struct, msh_side, href(x{:}));
    
        M2(1:ndofs,1:ndofs) = M2(1:ndofs,1:ndofs) + CC.' * op_gradu_n_gradv_n (sp_bnd_struct, sp_bnd_struct, msh_side, coeff_at_qnodes) * CC;
        rhs2(1:ndofs) = rhs2(1:ndofs) + CC.' * op_gradv_n_f (sp_bnd_struct, msh_side, dudn_at_qnodes); % I am missing the other part of the vector. It is in M2 :-)
      end
    end
  end
end

u_drchlt = M(drchlt_dofs, drchlt_dofs) \ rhs(drchlt_dofs, 1);

uu = sparse (hspace.ndof, 1);
uu(drchlt_dofs) = u_drchlt;

drchlt_dofs2 = setdiff (drchlt_dofs2, drchlt_dofs);
rhs2 = rhs2 - M2 * uu;
u_drchlt2 = M2(drchlt_dofs2, drchlt_dofs2) \ rhs2(drchlt_dofs2);

uu(drchlt_dofs2) = u_drchlt2;

drchlt_dofs = union (drchlt_dofs, drchlt_dofs2);
u_drchlt = uu(drchlt_dofs);

end