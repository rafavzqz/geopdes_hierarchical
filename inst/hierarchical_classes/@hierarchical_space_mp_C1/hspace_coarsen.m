% HSPACE_COARSEN: coarsen the hierarchical space, updating the fields of the object.
%
%   hspace = hspace_coarsen (hspace, hmsh, funs_to_reactivate, removed_cells)
%
% INPUT:
%
%   hspace:         object representing the fine hierarchical space (see hierarchical_space_mp_C1)
%   hmsh:           object representing the hierarchical mesh, already coarsened (see hierarchical_mesh_mp)
%   funs_to_reactivate: cell array with the indices of the functions to reactivate for each level
%   removed_cells:  cell array with the elements removed during coarsening, for each level
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

function hspace = hspace_coarsen (hspace, hmsh, FtR, removed_cells)

% Update active functions
hspace = update_active_functions (hspace, hmsh, FtR, removed_cells);

% Update the matrices for changing basis
[hspace.Csub, hspace.Csub_row_indices] = hspace_subdivision_matrix (hspace, hmsh);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function hspace = update_active_functions (hspace, hmsh, marked_funs, removed_cells)
%
% This function updates the active (hspace.active) and deactivated (hspace.deactivated) degrees of freedom,
% reactivating the functions in marked_funs. 
% The function also updates hspace.nlevels, hspace.ndof and hspace.ndof_per_level
%
% Input:    hspace:    the fine space, an object of the class hierarchical_space
%           hmsh:      an object of the class hierarchical_mesh, already coarsened
%           marked_funs{lev}: indices of active functions of level lev to be reactivated
%           removed_cells{lev}: indices of cells of level lev that have been removed 
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

% if (hspace.truncated)
%   hspace.coeff_pou = ones (hspace.ndof, 1);
% else
%   hspace.coeff_pou = Cref * hspace.coeff_pou;
% end

end
