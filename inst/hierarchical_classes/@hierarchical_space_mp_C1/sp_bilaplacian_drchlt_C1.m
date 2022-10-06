function [u_drchlt, drchlt_dofs, kernel_info] = sp_bilaplacian_drchlt_C1 (hspace, hmsh, refs, h, dudn)

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
        
% Since "charlen" is not present in msh_side, I assume isotropic elements
        charlen = msh_bnd_struct.element_size;

%        M(1:ndofs,1:ndofs) = M(1:ndofs,1:ndofs) + CC(sp_bnd.dofs,:).' * op_u_v_tp (sp_bnd, sp_bnd, msh_bnd, coeff_at_qnodes) * CC(sp_bnd.dofs,:);
        M(1:ndofs,1:ndofs) = M(1:ndofs,1:ndofs) + CC.' * op_u_v (sp_bnd_struct, sp_bnd_struct, msh_side, coeff_at_qnodes) * CC;
        rhs(1:ndofs) = rhs(1:ndofs) + CC.' * op_f_v (sp_bnd_struct, msh_side, href(x{:}));
    
        M2(1:ndofs,1:ndofs) = M2(1:ndofs,1:ndofs) + CC.' * op_gradu_n_gradv_n (sp_bnd_struct, sp_bnd_struct, msh_side, coeff_at_qnodes.*charlen) * CC;
        rhs2(1:ndofs) = rhs2(1:ndofs) + CC.' * op_gradv_n_f (sp_bnd_struct, msh_side, dudn_at_qnodes.*charlen); % I am missing the other part of the vector. It is in M2 :-)
      end
    end
  end
end

M_bdry = M + M2;
dofs_to_remove = [];
vertices_numbers = [];
row_indices = [];
count_vert = 0;

% Check the kernel of vertex functions on Dirichlet boundary vertices
% Pick up the basis function with the max. abs. coeff in the kernel, 
%  remove it from drchlt_dofs, and add the function in the kernel into the
%  internal part (it goes in the output)
B_change_local = [];
n_boundaries = numel(hmsh.mesh_of_level(1).boundaries); % number of boundary edges
global_refs = numel(hspace.space_of_level(1).interfaces) - n_boundaries + refs; % global numbering of Dirichlet boundary edges

for iv = 1 : numel(hspace.space_of_level(1).vertices) % Loop on the vertices
    % Loop just over Dirichlet boundary vertices
    if ~isempty(intersect(global_refs, hspace.space_of_level(1).vertices(iv).edges))
        if (hspace.space_of_level(1).vertices(iv).boundary_vertex)
            for lev = 1 : hspace.nlevels % Loop on the levels    
                [~, ind_active_vert] = intersect(hspace.active{lev}, hspace.space_of_level(lev).dofs_on_vertex{iv}); % Indices of the functions of vertex iv in the vector of active functions of level lev
                if ~isempty(ind_active_vert) % Check if the vertex iv is active on level lev
                    if numel(hspace.space_of_level(lev).vertices(iv).patches) > 1 % Check if the vertex is shared among more than one patch
                        
                        patches = hspace.space_of_level(lev).vertices(iv).patches([1 end]);

                        operations = hspace.space_of_level(lev).vertices(iv).patch_reorientation([1 end], :);
                        indices_loc_R = indices_reorientation(hspace.space_of_level(lev).sp_patch{patches(1)}.ndof_dir, operations(1, :));
                        indices_loc_L = indices_reorientation(hspace.space_of_level(lev).sp_patch{patches(2)}.ndof_dir, operations(2, :));

                        indices_loc_R = indices_loc_R(:);
                        indices_loc_L = indices_loc_L(:);

                        Cpatch_ind_R = indices_loc_R([2 3 hspace.space_of_level(lev).sp_patch{patches(1)}.ndof_dir(1)+2]);
                        Cpatch_ind_L = indices_loc_L([hspace.space_of_level(lev).sp_patch{patches(1)}.ndof_dir(1)+1 hspace.space_of_level(lev).sp_patch{patches(1)}.ndof_dir(1)+2 2*hspace.space_of_level(lev).sp_patch{patches(1)}.ndof_dir(1)+1]);

                        M_ker = [hspace.space_of_level(lev).Cpatch{patches(1)}(Cpatch_ind_R, hspace.space_of_level(lev).dofs_on_vertex{iv}); ...
                                 hspace.space_of_level(lev).Cpatch{patches(2)}(Cpatch_ind_L, hspace.space_of_level(lev).dofs_on_vertex{iv})];

                        ker = null(full(M_ker));
                            if (~isempty(ker))
                                count_vert = count_vert + 1;
                                [~, ind] = max(abs(ker));

                                row_inds = (count_vert-1)*6 + (1:6);
                                B_change_local = blkdiag (B_change_local, ker(:));

                                dofs_on_vertex = sum(hspace.ndof_per_level(1:lev-1)) + ind_active_vert;
                                vertices_numbers(count_vert) = iv;
                                dofs_to_remove(count_vert) = dofs_on_vertex(ind);
                                row_indices(row_inds) = dofs_on_vertex;
                            end
                    end
                end
            end
        end
    end
end
kernel_info = struct ('vertices_numbers', vertices_numbers, 'all_vertex_dofs', row_indices, 'quasi_interior_dofs', dofs_to_remove, 'B_change_local', sparse(B_change_local));

drchlt_dofs = union (drchlt_dofs, drchlt_dofs2);
drchlt_dofs = setdiff(drchlt_dofs, dofs_to_remove);

u_drchlt = M_bdry(drchlt_dofs,drchlt_dofs) \ (rhs(drchlt_dofs) + rhs2(drchlt_dofs));

end

function indices = indices_reorientation (ndof_dir, operations)
  ndof = prod (ndof_dir);
  indices = reshape (1:ndof, ndof_dir);
  if (operations(1))
    indices = flipud (indices);
  end
  if (operations(2))
    indices = fliplr (indices);
  end
  if (operations(3))
    indices = indices.';
  end   
end


% u_drchlt = M(drchlt_dofs, drchlt_dofs) \ rhs(drchlt_dofs, 1);
% 
% uu = sparse (hspace.ndof, 1);
% uu(drchlt_dofs) = u_drchlt;
% 
% drchlt_dofs2 = setdiff (drchlt_dofs2, drchlt_dofs);
% rhs2 = rhs2 - M2 * uu;
% u_drchlt2 = M2(drchlt_dofs2, drchlt_dofs2) \ rhs2(drchlt_dofs2);
% 
% uu(drchlt_dofs2) = u_drchlt2;
% 
% drchlt_dofs = union (drchlt_dofs, drchlt_dofs2);
% u_drchlt = uu(drchlt_dofs);
% 
% end