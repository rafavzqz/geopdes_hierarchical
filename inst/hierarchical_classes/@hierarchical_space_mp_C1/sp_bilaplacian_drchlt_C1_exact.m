% SP_BILAPLACIAN_DRCHLT_C1: assign the degrees of freedom of essential boundary conditions (value and normal derivative) through a projection.
%  On boundary vertices, the kernel is computed to remove linear dependencies when restricting the functions to the boundary.
%
%   [u, dofs, kernel_info] = sp_bilaplacian_drchlt_C1 (hspace, hmsh, refs, h, dudn)
%
% INPUT:
%
%  hspace:     object representing the space of hierarchical splines (see hierarchical_space_mp_C1)
%  hmsh:       object representing the hierarchical mesh (see hierarchical_mesh_mp)
%  refs:       boundary references on which the conditions are imposed
%  h:          function handle to compute the Dirichlet condition
%  dudn:       function handle to compute the Neumann condition
%
% OUTPUT:
%
%  u:           assigned value to the degrees of freedom
%  dofs:        global numbering of the corresponding basis functions
%  kernel_info: a struct with information of kernel computation, containing:
%              - vertices_numbers: vertices which contain a function in the kernel
%              - all_vertex_dofs:  all functions on those vertices
%              - quasi_interior_dofs: functions that will be treated as
%                                     internal ones (as many as in the kernel)
%              - B_change_local: coefficients of the functions in the kernel,
%                                in terms of vertex basis functions. Matrix of size
%                                numel(all_vertex_dofs) x numel (quasi_interior_dofs)
%
% Copyright (C) 2022-2023 Cesare Bracco, Andrea Farahat, Rafael Vazquez
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

function [u_drchlt, drchlt_dofs, kernel_info] = sp_bilaplacian_drchlt_C1_exact (hspace, hmsh, refs, uex, gradex)

M = spalloc (hspace.ndof, hspace.ndof, hspace.ndof);
rhs = zeros (hspace.ndof, 1);

M2 = spalloc (hspace.ndof, hspace.ndof, hspace.ndof);
rhs2 = zeros (hspace.ndof, 1);

drchlt_dofs = [];
drchlt_dofs2 = [];

boundaries = hmsh.mesh_of_level(1).boundaries;
for iref = refs
%   href = @(varargin) h(varargin{:}, iref);
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
        msh_side = msh_eval_boundary_side (msh_lev_patch, iside, active_cells_patch);
        msh_side_from_interior = msh_boundary_side_from_interior (msh_lev_patch, iside);

        sp_bnd = sp_lev_patch.constructor (msh_side_from_interior);
        msh_bnd_struct = msh_evaluate_element_list (msh_side_from_interior, active_cells_patch);
        sp_bnd_struct = sp_evaluate_element_list (sp_bnd, msh_bnd_struct, 'value', true, 'gradient', true);

        sp_bnd = sp_lev_patch.boundary(iside);
        [Cpatch, Cpatch_cols_lev] = sp_compute_Cpatch (hspace.space_of_level(ilev), iptc);
        [~,Csub_rows,Cpatch_cols] = intersect (hspace.Csub_row_indices{ilev}, Cpatch_cols_lev);
        CC = Cpatch(:,Cpatch_cols) * hspace.Csub{ilev}(Csub_rows,:);

        [~,icol] = find (CC(sp_bnd.dofs,:));
        [~,jcol] = find (CC(sp_bnd.adjacent_dofs,:));
    
        drchlt_dofs = union (drchlt_dofs, icol);
        drchlt_dofs2 = union (drchlt_dofs2, jcol);
    
        for idim = 1:hmsh.rdim
          x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
        end
        coeff_at_qnodes = ones (size(x{1}));
        dudn_at_qnodes = reshape (sum (gradex(x{:}) .* msh_side.normal, 1), msh_side.nqn, msh_side.nel);
        
% Since "charlen" is not present in msh_side, I assume isotropic elements
        charlen = msh_bnd_struct.element_size;

%        M(1:ndofs,1:ndofs) = M(1:ndofs,1:ndofs) + CC(sp_bnd.dofs,:).' * op_u_v_tp (sp_bnd, sp_bnd, msh_bnd, coeff_at_qnodes) * CC(sp_bnd.dofs,:);
        M(1:ndofs,1:ndofs) = M(1:ndofs,1:ndofs) + CC.' * op_u_v (sp_bnd_struct, sp_bnd_struct, msh_side, coeff_at_qnodes) * CC;
        rhs(1:ndofs) = rhs(1:ndofs) + CC.' * op_f_v (sp_bnd_struct, msh_side, uex(x{:}));
    
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
  if (~isempty(intersect(global_refs, hspace.space_of_level(1).vertices(iv).edges)))
    if (hspace.space_of_level(1).vertices(iv).boundary_vertex)
      for lev = 1 : hspace.nlevels % Loop on the levels    
        [~, ind_active_vert] = intersect(hspace.active{lev}, hspace.space_of_level(lev).dofs_on_vertex{iv}); % Indices of the functions of vertex iv in the vector of active functions of level lev
        if (~isempty(ind_active_vert)) % Check if the vertex iv is active on level lev
          if (numel(hspace.space_of_level(lev).vertices(iv).patches) > 1) % Check if the vertex is shared among more than one patch
                        
             patches = hspace.space_of_level(lev).vertices(iv).patches([1 end]);

             operations = hspace.space_of_level(lev).vertices(iv).patch_reorientation([1 end], :);
             indices_loc_R = indices_reorientation(hspace.space_of_level(lev).sp_patch{patches(1)}.ndof_dir, operations(1, :));
             indices_loc_L = indices_reorientation(hspace.space_of_level(lev).sp_patch{patches(2)}.ndof_dir, operations(2, :));

             indices_loc_R = indices_loc_R(:);
             indices_loc_L = indices_loc_L(:);

             Cpatch_ind_R = indices_loc_R([2 3 hspace.space_of_level(lev).sp_patch{patches(1)}.ndof_dir(1)+2]);
             Cpatch_ind_L = indices_loc_L([hspace.space_of_level(lev).sp_patch{patches(1)}.ndof_dir(1)+1 hspace.space_of_level(lev).sp_patch{patches(1)}.ndof_dir(1)+2 2*hspace.space_of_level(lev).sp_patch{patches(1)}.ndof_dir(1)+1]);

            [Cpatch1, Cpatch_cols1] = sp_compute_Cpatch (hspace.space_of_level(lev), patches(1));
            [Cpatch2, Cpatch_cols2] = sp_compute_Cpatch (hspace.space_of_level(lev), patches(2));
            [~,~,inds1] = intersect (hspace.space_of_level(lev).dofs_on_vertex{iv}, Cpatch_cols1);
            [~,~,inds2] = intersect (hspace.space_of_level(lev).dofs_on_vertex{iv}, Cpatch_cols2);

            M_ker = [Cpatch1(Cpatch_ind_R, inds1); ...
                     Cpatch2(Cpatch_ind_L, inds2)];

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
