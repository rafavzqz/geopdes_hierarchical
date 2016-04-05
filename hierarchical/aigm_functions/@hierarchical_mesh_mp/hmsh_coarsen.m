% HMSH_COARSEN: coarsen the hierarchical mesh, updating the fields of the object.
%
%   [hmsh, new_elements] = hmsh_coarsen (hmsh, marked)
%
% INPUT:
%
%   hmsh:    object representing the fine hierarchical mesh (see hierarchical_mesh_mp)
%   marked:  cell array with the indices, in the Cartesian grid, of the elements to be reactivated for each level
%
% OUTPUT:
%
%   hmsh:             object representing the coarsened hierarchical mesh (see hierarchical_mesh_mp)
%   removed_elements: cell array with the indices of the elements removed for each level
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
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

function [hmsh, removed_elements] = hmsh_coarsen (hmsh, M)

boundary = ~isempty (hmsh.boundary);

% Update the set of active elements
old_elements = hmsh.active;
[hmsh, removed_elements] = update_active_cells (hmsh, M);

% Update msh_lev
hmsh.msh_lev = update_msh_lev (hmsh, old_elements);
% hmsh.msh_lev = cell (hmsh.nlevels,1);
% for ilev = 1 : hmsh.nlevels
%   hmsh.msh_lev{ilev} = msh_evaluate_element_list (hmsh.mesh_of_level(ilev), hmsh.active{ilev});
% end


% Update the boundary, calling the function recursively
if (boundary)
  if (hmsh.ndim > 1)
    M_boundary = cell (size (M));
    for lev = 1:min(numel (M), hmsh.boundary.nlevels)
      Nelem = cumsum ([0 hmsh.mesh_of_level(lev).nel_per_patch]);
      Nelem_bnd = cumsum ([0 hmsh.boundary.mesh_of_level(lev).nel_per_patch]);
      for iptc = 1:hmsh.boundary.npatch
        patch_number = hmsh.boundary.mesh_of_level(1).patch_numbers(iptc);
        side_number  = hmsh.boundary.mesh_of_level(1).side_numbers(iptc);
        [~,indices,~] = intersect (Nelem(patch_number)+1:Nelem(patch_number+1), M{lev});
        bnd_indices = get_boundary_indices (side_number, hmsh.mesh_of_level(lev).msh_patch{patch_number}.nel_dir, indices);
        M_boundary{lev} = union (M_boundary{lev}, bnd_indices+Nelem_bnd(iptc));
      end
    end
    hmsh.boundary = hmsh_coarsen (hmsh.boundary, M_boundary);
  end
else
  hmsh.boundary = [];
end

% Remove last level if empty (this is not a requirement)
if (isempty (hmsh.active{end}))
  hmsh = remove_empty_level (hmsh);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hmsh, removed_cells] = update_active_cells (hmsh, M)
%
% function [hmsh, removed_cells] = update_active_cells (hmsh, M)
%
% Update the sets of active cells (hmsh.active) and deactivated cells (hmsh.deactivated) in each level when
% reactivating the cells in M. This function also updates hmsh.nlevels, hmsh.nel and hmsh.nel_per_level
%
% INPUT
%     hmsh: object representing the fine hierarchical mesh (see hierarchical_mesh_mp)
%     M{lev}: global indices of marked cells of level lev (one row per cell)
%
% OUTPUT
%     hmsh: object representing the coarsened hierarchical mesh (see hierarchical_mesh_mp)
%     removed_cells{lev}: global indices of the removed cells of level lev (one row per cell)
%

nlevels = hmsh.nlevels;
removed_cells = cell (nlevels+1, 1);

for lev = nlevels-1:-1:1
  removed_cells{lev+1} = hmsh_get_children (hmsh, lev, M{lev});
  hmsh.active{lev+1} = setdiff (hmsh.active{lev+1}, removed_cells{lev+1});
  hmsh.deactivated{lev} = setdiff (hmsh.deactivated{lev}, M{lev});
  hmsh.active{lev} = union (hmsh.active{lev}, M{lev});
end

if (isempty (hmsh.active{nlevels}))
  hmsh.nlevels = nlevels-1;
  hmsh.active = hmsh.active(1:nlevels-1);
  hmsh.deactivated = hmsh.deactivated(1:nlevels-1);
end

% Update hmsh.nel_per_level and hmsh.nel
hmsh.nel_per_level = cellfun (@numel, hmsh.active);
hmsh.nel = sum (hmsh.nel_per_level);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msh_lev = update_msh_lev (hmsh, old_elements)
%
% function msh_lev = update_msh_lev (hmsh, old_elements)
%
% Update the information in msh_lev, computing only the elements that have been added to the mesh
%
% INPUT
%     hmsh: object representing the fine hierarchical mesh (see hierarchical_mesh_mp)
%     old_elements{lev}:     active elements in the previous mesh
%
% OUTPUT
%     msh_lev: the structures with the msh information for each level
%

msh_lev = cell (hmsh.nlevels, 1);

for lev = 1:hmsh.nlevels
  if (numel (old_elements{lev}) == 0)
    msh_lev{lev} = msh_evaluate_element_list (hmsh.mesh_of_level(lev), hmsh.active{lev});
  else
    new_elements = setdiff (hmsh.active{lev}, old_elements{lev});
      
    [~, iold, iold_act] = intersect (old_elements{lev}, hmsh.active{lev});
    msh_lev{lev}.npatch = hmsh.npatch;
    msh_lev{lev}.ndim = hmsh.ndim;
    msh_lev{lev}.rdim = hmsh.rdim;
    msh_lev{lev}.nel = hmsh.nel_per_level(lev);
    msh_lev{lev}.elem_list = hmsh.active{lev}(:).';

    if (isempty (new_elements))
      indices = iold_act;
      msh_new = struct ('quad_weights', [], 'geo_map', [], 'geo_map_jac', [], 'geo_map_der2', [], 'jacdet', [], 'element_size', []);
    else
      msh_new = msh_evaluate_element_list (hmsh.mesh_of_level(lev), new_elements);
      [~, ~, inew_act] = intersect (new_elements, hmsh.active{lev});
      indices = [iold_act(:); inew_act(:)];
    end
    msh_lev{lev}.quad_weights(:,indices) = [hmsh.msh_lev{lev}.quad_weights(:,iold), msh_new.quad_weights];
    msh_lev{lev}.geo_map(:,:,indices) = cat (3, hmsh.msh_lev{lev}.geo_map(:,:,iold), msh_new.geo_map);
    msh_lev{lev}.geo_map_jac(:,:,:,indices) = cat (4, hmsh.msh_lev{lev}.geo_map_jac(:,:,:,iold), msh_new.geo_map_jac);
    msh_lev{lev}.geo_map_der2(:,:,:,:,indices) = cat (5, hmsh.msh_lev{lev}.geo_map_der2(:,:,:,:,iold), msh_new.geo_map_der2);
    msh_lev{lev}.jacdet(:,indices) = [hmsh.msh_lev{lev}.jacdet(:,iold), msh_new.jacdet];
    msh_lev{lev}.element_size(:,indices) = [hmsh.msh_lev{lev}.element_size(:,iold), msh_new.element_size];

    Nelem = cumsum ([0 hmsh.mesh_of_level(lev).nel_per_patch]);
    for iptc = 1:msh_lev{lev}.npatch
      [~,indices,~] = intersect (Nelem(iptc)+1:Nelem(iptc+1), msh_lev{lev}.elem_list);
      msh_lev{lev}.nel_per_patch(iptc) = numel (indices);
      msh_lev{lev}.elem_list_of_patch{iptc} = indices;
      msh_lev{lev}.nel_dir_of_patch{iptc} = hmsh.mesh_of_level(lev).msh_patch{iptc}.nel_dir;
    end
    msh_lev{lev}.nqn = hmsh.mesh_of_level(lev).msh_patch{1}.nqn;
    msh_lev{lev}.nqn_dir = hmsh.mesh_of_level(lev).msh_patch{1}.nqn_dir;
  end
end

end
