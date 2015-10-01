% HMSH_REFINE: refine the hierarchical mesh, updating the fields of the object.
%
%   [hmsh, new_elements] = hmsh_refine (hmsh, marked, indices)
%
% INPUT:
%
%   hmsh:    object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   marked:  cell array with the indices, in the Cartesian grid, of the marked elements for each level
%   indices: relative position of the marked elements in the numbering of the hierarchical mesh
%         within the level, that is, in hmsh.active{lev}
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

function [hmsh, new_elements] = hmsh_refine (hmsh, M, indices)

boundary = ~isempty (hmsh.boundary);

% Computation of a Cartesian grid if a new level is activated
if (~isempty(M{hmsh.nlevels}))
  rule = msh_gauss_nodes (hmsh.mesh_of_level(1).nqn_dir);
  [dummy, zeta] = kntrefine (hmsh.mesh_of_level(hmsh.nlevels).breaks, hmsh.nsub-1, ones(1, hmsh.ndim), zeros(1, hmsh.ndim));
  [qn, qw] = msh_set_quad_nodes (zeta, rule);
    
  hmsh.mesh_of_level(hmsh.nlevels+1) = msh_cartesian (zeta, qn, qw, hmsh.geometry, 'boundary', boundary);
end

% Update the set of active elements
[hmsh, new_elements] = update_active_cells (hmsh, M, indices);

% Update msh_lev
% XXXX This can be improved to save computational time, avoiding to
% recompute the elements that did not change
hmsh.msh_lev = cell (hmsh.nlevels,1);
for ilev = 1 : hmsh.nlevels
  hmsh.msh_lev{ilev} = msh_evaluate_element_list (hmsh.mesh_of_level(ilev), hmsh.active{ilev});
end


% Update the boundary , calling the function recursively
if (boundary)
  if (hmsh.ndim > 1)
    for iside = 1:2*hmsh.ndim
      M_boundary = cell (size (M));
      for lev = 1:numel (M)
        M_boundary{lev} = get_boundary_indices (iside, hmsh.mesh_of_level(lev).nel_dir, M{lev});
      end
      hmsh.boundary(iside) = hmsh_refine (hmsh.boundary(iside), M_boundary, []);
    end
  end
else
  hmsh.boundary = [];
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hmsh, new_cells] = update_active_cells (hmsh, M, indices)
%
% function [hmsh, new_cells] = update_active_cells (hmsh, M, indices)
%
% This function updates the active cells (hmsh.active and hmsh.globnum_active) and deactivated cells (hmsh.deactivated) in each level when
% refining the cells in M. This function also updates hmsh.nlevels, hmsh.nel and hmsh.nel_per_level
%
% INPUT
%     hmsh: object representing the coarse hierarchical mesh (see hierarchical_mesh)
%     M{lev}: global indices of marked cells of level lev (one row per cell)
%     indices{lev}: array such that M{lev} =
%           hmsh.active{lev}(indices{lev}). If indices is [], then this function computes this information.
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
    if (isempty(indices))
      [dummy, indE] = ismember (M{lev}, hmsh.active{lev});
%       if (~all (dummy))
%         warning('update_active_cells: Some nonactive cells were selected');
%       end
    else
      indE = indices{lev};
    end
    hmsh.active{lev}(indE) = [];
      
    new_cells{lev+1} = split_cells_of_level (hmsh, lev, M{lev});
    hmsh.deactivated{lev} = union (hmsh.deactivated{lev}, M{lev});
  end
end

% Update the active cells, by adding the children of the refined cells
for lev = 1:nlevels
  if (~isempty(M{lev}))
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
  children = union (children, sub2ind ([hmsh.mesh_of_level(lev+1).nel_dir, 1], z{:}));
end

end
