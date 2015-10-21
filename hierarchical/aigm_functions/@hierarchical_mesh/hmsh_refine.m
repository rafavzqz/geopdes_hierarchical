% HMSH_REFINE: refine the hierarchical mesh, updating the fields of the object.
%
%   [hmsh, new_elements] = hmsh_refine (hmsh, marked)
%
% INPUT:
%
%   hmsh:    object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   marked:  cell array with the indices, in the Cartesian grid, of the marked elements for each level
%
% OUTPUT:
%
%   hmsh:         object representing the refined hierarchical mesh (see hierarchical_mesh)
%   new_elements: cell array with the global indices of the new active elements for each level
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
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

function [hmsh, new_elements] = hmsh_refine (hmsh, M)

boundary = ~isempty (hmsh.boundary);

% Computation of a Cartesian grid if a new level is activated
if (~isempty(M{hmsh.nlevels}))
  rule = msh_gauss_nodes (hmsh.mesh_of_level(1).nqn_dir);
  [dummy, zeta] = kntrefine (hmsh.mesh_of_level(hmsh.nlevels).breaks, hmsh.nsub-1, ones(1, hmsh.ndim), zeros(1, hmsh.ndim));
  [qn, qw] = msh_set_quad_nodes (zeta, rule);
    
  hmsh.mesh_of_level(hmsh.nlevels+1) = msh_cartesian (zeta, qn, qw, hmsh.geometry, 'boundary', boundary);
  hmsh.mesh_of_level(hmsh.nlevels+1).boundary = []; % Remove redundant information
end

% Update the set of active elements
old_elements = hmsh.active;
[hmsh, new_elements] = update_active_cells (hmsh, M);

% Update msh_lev
hmsh.msh_lev = update_msh_lev (hmsh, old_elements, new_elements);
% hmsh.msh_lev = cell (hmsh.nlevels,1);
% for ilev = 1 : hmsh.nlevels
%   hmsh.msh_lev{ilev} = msh_evaluate_element_list (hmsh.mesh_of_level(ilev), hmsh.active{ilev});
% end

% Update the boundary , calling the function recursively
if (boundary)
  if (hmsh.ndim > 1)
    for iside = 1:2*hmsh.ndim
      M_boundary = cell (size (M));
      for lev = 1:numel (M)
        M_boundary{lev} = get_boundary_indices (iside, hmsh.mesh_of_level(lev).nel_dir, M{lev});
      end
      hmsh.boundary(iside) = hmsh_refine (hmsh.boundary(iside), M_boundary);
    end
  end
else
  hmsh.boundary = [];
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hmsh, new_cells] = update_active_cells (hmsh, M)
%
% function [hmsh, new_cells] = update_active_cells (hmsh, M)
%
% This function updates the active cells (hmsh.active and hmsh.globnum_active) and deactivated cells (hmsh.deactivated) in each level when
% refining the cells in M. This function also updates hmsh.nlevels, hmsh.nel and hmsh.nel_per_level
%
% INPUT
%     hmsh: object representing the coarse hierarchical mesh (see hierarchical_mesh)
%     M{lev}: global indices of marked cells of level lev (one row per cell)
%
% OUTPUT
%     hmsh: object representing the refined hierarchical mesh (see hierarchical_mesh)
%     new_cells{lev}: global indices of the new cells of level lev (one row per cell)
%

if (nargin == 2)
  indices = [];
end

nlevels = hmsh.nlevels;

if (~isempty(M{nlevels})) % if a new level is going to be activated
  hmsh.nlevels = hmsh.nlevels + 1;
  hmsh.active{nlevels+1} = [];
  hmsh.deactivated{nlevels+1} = [];
  new_cells = cell (nlevels+1, 1);
else
  new_cells = cell (nlevels, 1);
end

% Deactivate the cells to be refined, and compute their children
for lev = 1:nlevels
  if (~isempty(M{lev}))
    [dummy, indE] = ismember (M{lev}, hmsh.active{lev});
%       if (~all (dummy))
%         warning('update_active_cells: Some nonactive cells were selected');
%       end
    hmsh.active{lev}(indE) = [];
    hmsh.deactivated{lev} = union (hmsh.deactivated{lev}, M{lev});
      
    new_cells{lev+1} = split_cells_of_level (hmsh, lev, M{lev});
    hmsh.active{lev+1} = union (hmsh.active{lev+1}, new_cells{lev+1});
  end
end

% Update hmsh.nel_per_level and hmsh.nel
hmsh.nel_per_level = cellfun (@numel, hmsh.active);
hmsh.nel = sum (hmsh.nel_per_level);

% Update hmsh.globnum_active
hmsh.globnum_active = zeros (hmsh.nel,hmsh.ndim+1);
Ne = cumsum ([0; hmsh.nel_per_level(:)]);
for lev = 1:hmsh.nlevels
  ind_e = (Ne(lev)+1):Ne(lev+1);
  if (~isempty(ind_e))
    hmsh.globnum_active(ind_e, 1) = lev;
    globnum = cell (1,hmsh.ndim);
    [globnum{:}] = ind2sub (hmsh.mesh_of_level(lev).nel_dir, hmsh.active{lev}(:));
    hmsh.globnum_active(ind_e, 2:end) = cell2mat (globnum);
  end
end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function children = split_cells_of_level (hmsh, lev, ind)
%
% function children = split_cells_of_level (hmsh, lev, ind)
%
% Split a set of cells of hmsh, of level lev, with the subdivision given by hmsh.nsub.
%
% Input:
%     hmsh: the hierarchical mesh
%     lev:  level of the cells to subdivide
%     ind:  indices of the cells in the Cartesian grid
%
% Output:
%     children: indices of the children, with the numbering of the Cartesian grid
%

z = cell (hmsh.ndim, 1);
cells_sub = cell (hmsh.ndim, 1);
[cells_sub{:}] = ind2sub ([hmsh.mesh_of_level(lev).nel_dir, 1], ind); % The extra 1 makes it work in any dimension

children = [];
for ii = 1:numel(cells_sub{1})
  aux = cell (hmsh.ndim, 1);
  for idim = 1:hmsh.ndim
    aux{idim} = hmsh.nsub(idim)*(cells_sub{idim}(ii)-1)+1:hmsh.nsub(idim)*(cells_sub{idim}(ii));
  end
  [z{1:hmsh.ndim}] = ndgrid (aux{:});
  auxI = sub2ind ([hmsh.mesh_of_level(lev+1).nel_dir, 1], z{:});
  children = union (children, auxI(:));
end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msh_lev = update_msh_lev (hmsh, old_elements, new_elements)
%
% function msh_lev = update_msh_lev (hmsh, old_elements, new_elements)
%
% This function updates the information in msh_lev, computing only the elements that have been added to the mesh
%
% INPUT
%     hmsh: object representing the fine hierarchical mesh (see hierarchical_mesh)
%     old_elements{lev}: active elements in the previous coarse mesh
%     new_elements{lev}: elements that have been added after refinement
%
% OUTPUT
%     msh_lev: the structures with the msh information for each level
%

msh_lev = cell (hmsh.nlevels, 1);

for lev = 1:hmsh.nlevels
  if (lev > numel (old_elements) || numel (old_elements{lev}) == 0)
    msh_lev{lev} = msh_evaluate_element_list (hmsh.mesh_of_level(lev), hmsh.active{lev});
  else
    [~, iold, iold_act] = intersect (old_elements{lev}, hmsh.active{lev});
    msh_lev{lev}.ndim = hmsh.ndim;
    msh_lev{lev}.rdim = hmsh.rdim;
    msh_lev{lev}.nel = hmsh.nel_per_level(lev);
    msh_lev{lev}.elem_list = hmsh.active{lev}(:)';
    msh_lev{lev}.nel_dir = hmsh.mesh_of_level(lev).nel_dir;
    msh_lev{lev}.nqn_dir = hmsh.mesh_of_level(lev).nqn_dir;
    msh_lev{lev}.nqn = hmsh.mesh_of_level(lev).nqn;

    if (isempty (new_elements{lev}))
      indices = iold_act;
      msh_new = struct ('quad_weights', [], 'geo_map', [], 'geo_map_jac', [], 'geo_map_der2', [], 'jacdet', [], 'element_size', []);
    else
      msh_new = msh_evaluate_element_list (hmsh.mesh_of_level(lev), new_elements{lev});
      [~, ~, inew_act] = intersect (new_elements{lev}, hmsh.active{lev});
      indices = [iold_act(:); inew_act(:)];
    end
    msh_lev{lev}.quad_weights(:,indices) = [hmsh.msh_lev{lev}.quad_weights(:,iold), msh_new.quad_weights];
    msh_lev{lev}.geo_map(:,:,indices) = cat (3, hmsh.msh_lev{lev}.geo_map(:,:,iold), msh_new.geo_map);
    msh_lev{lev}.geo_map_jac(:,:,:,indices) = cat (4, hmsh.msh_lev{lev}.geo_map_jac(:,:,:,iold), msh_new.geo_map_jac);
    msh_lev{lev}.geo_map_der2(:,:,:,:,indices) = cat (5, hmsh.msh_lev{lev}.geo_map_der2(:,:,:,:,iold), msh_new.geo_map_der2);
    msh_lev{lev}.jacdet(:,indices) = [hmsh.msh_lev{lev}.jacdet(:,iold), msh_new.jacdet];
    msh_lev{lev}.element_size(:,indices) = [hmsh.msh_lev{lev}.element_size(:,iold), msh_new.element_size];
  end
end

end
