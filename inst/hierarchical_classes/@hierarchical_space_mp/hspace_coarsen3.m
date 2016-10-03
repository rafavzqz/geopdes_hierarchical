% HSPACE_COARSEN: coarsen the hierarchical space, updating the fields of the object.
%
%   hspace = hspace_coarsen (hspace, hmsh, functions_to_remove)
%
% INPUT:
%
%   hspace:         object representing the fine hierarchical space (see hierarchical_space_mp)
%   hmsh:           object representing the hierarchical mesh, already coarsened (see hierarchical_mesh_mp)
%   funs_to_remove: cell array with the indices, in the tensor product setting, of the functions to remove for each level
%
% OUTPUT:
%
%   hspace:    object representing the coarsened hierarchical space (see hierarchical_space_mp)
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

function hspace = hspace_coarsen3 (hspace, hmsh, FtR, removed_cells)

boundary = ~isempty (hspace.boundary);

% Update active functions
hspace = update_active_functions (hspace, hmsh, FtR, removed_cells);

% Update the matrices for changing basis
hspace.Csub = hspace_subdivision_matrix (hspace, hmsh);

% Fill the information for the boundaries
if (boundary)% && hmsh.ndim > 1)
  if (hmsh.ndim > 1)
    FtR_boundary = cell (size (FtR));
    levels = find (~cellfun (@isempty, FtR));
    for lev = levels(:).'
      [~,~,FtR_boundary{lev}] = intersect (FtR{lev}, hspace.space_of_level(lev).boundary.dofs);
    end
    cells_boundary = cell (size (removed_cells));
    for lev = 1:numel (removed_cells)
      Nelem = cumsum ([0 hmsh.mesh_of_level(lev).nel_per_patch]);
      Nelem_bnd = cumsum ([0 hmsh.boundary.mesh_of_level(lev).nel_per_patch]);
      for iptc = 1:hmsh.boundary.npatch
        patch_number = hmsh.boundary.mesh_of_level(1).patch_numbers(iptc);
        side_number  = hmsh.boundary.mesh_of_level(1).side_numbers(iptc);
        [~,indices,~] = intersect (Nelem(patch_number)+1:Nelem(patch_number+1), removed_cells{lev});
        bnd_indices = get_boundary_indices (side_number, hmsh.mesh_of_level(lev).msh_patch{patch_number}.nel_dir, indices);
        cells_boundary{lev} = union (cells_boundary{lev}, bnd_indices+Nelem_bnd(iptc));
      end
    end
    
    hspace.boundary = hspace_coarsen3 (hspace.boundary, hmsh.boundary, FtR_boundary, cells_boundary);

    Nf = cumsum ([0, hspace.ndof_per_level]);
    dofs = [];
    for lev = 1:hspace.boundary.nlevels
      [~,iact,jact] = intersect (hspace.active{lev}, hspace.space_of_level(lev).boundary.dofs);
      [~, reorder] = sort (jact);
      dofs = vertcat (dofs, Nf(lev) + iact(reorder));
    end
    hspace.boundary.dofs = dofs;
  elseif (hmsh.ndim == 1)
    error ('The 1D multipatch has not been implemented')
  end
else
    hspace.boundary = [];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function hspace = update_active_functions(hspace, hmsh, new_cells, marked_fun)
%
% This function updates the active dofs (hspace.active), their coefficients (hspace.coeff_pou) and deactivated dofs (hspace.deactivated) in each level when
% coarsening (removing) the functions in marked_fun. This function also updates hspace.ndof and hspace.ndof_per_level
%
% Input:    hspace:    the fine space, an object of the class hierarchical_space_mp
%           hmsh:      an object of the class hierarchical_mesh_mp, already coarsened
%           marked_fun{lev}: indices of active functions of level lev to be removed
%
% Output:   hspace:    the coarsened space, an object of the class hierarchical_space_mp
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

function hspace = update_active_functions (hspace, hmsh, funs_to_reactivate, removed_cells)

active = hspace.active;
deactivated = hspace.deactivated;

for lev = hspace.nlevels:-1:2
  if (strcmpi (hspace.type, 'standard') && ~isempty (removed_cells{lev}))
    removed_funs = sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), removed_cells{lev});
    active{lev} = setdiff (active{lev}, removed_funs);
    active{lev-1} = union (active{lev-1}, funs_to_reactivate{lev-1});
    deactivated{lev-1} = setdiff (deactivated{lev-1}, funs_to_reactivate{lev-1});
    
  elseif (strcmpi (hspace.type, 'simplified') && ~isempty (funs_to_reactivate{lev-1}))
    active{lev-1} = union (active{lev-1}, funs_to_reactivate{lev-1});
    deactivated{lev-1} = setdiff (deactivated{lev-1}, funs_to_reactivate{lev-1});
    children = hspace_get_children (hspace, lev-1, funs_to_reactivate{lev-1});

    neighbors = sp_get_neighbors (hspace.space_of_level(lev-1), hmsh.mesh_of_level(lev-1), funs_to_reactivate{lev-1});
    deact_neighs = intersect (deactivated{lev-1}, neighbors);
    children = setdiff (children, hspace_get_children (hspace, lev-1, deact_neighs));
    active{lev} = setdiff (active{lev}, children);
  end
end

hspace.active = active(1:hspace.nlevels);
hspace.deactivated = deactivated(1:hspace.nlevels);
hspace.ndof_per_level = cellfun (@numel, hspace.active);
hspace.ndof = sum (hspace.ndof_per_level);

if (hspace.truncated)
  hspace.coeff_pou = ones (hspace.ndof, 1);
% else
%   hspace.coeff_pou = Cref * hspace.coeff_pou;
end

end
