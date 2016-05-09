% HSPACE_COARSEN: coarsen the hierarchical space, updating the fields of the object.
%
%   hspace = hspace_coarsen (hspace, hmsh, functions_to_remove)
%
% INPUT:
%
%   hspace:         object representing the fine hierarchical space (see hierarchical_space)
%   hmsh:           object representing the hierarchical mesh, already coarsened (see hierarchical_mesh)
%   funs_to_reactivate: cell array with the indices, in the tensor product setting, of the functions to reactivate for each level
%   removed_cells:  cell array with the elements removed during coarsening, for each level
%
% OUTPUT:
%
%   hspace:    object representing the coarsened hierarchical space (see hierarchical_space)
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

function hspace = hspace_coarsen2 (hspace, hmsh, FtR, removed_cells)

boundary = ~isempty (hspace.boundary);

% Update active functions
hspace = update_active_functions (hspace, hmsh, FtR, removed_cells);

% Update the matrices for changing basis
hspace.Csub = hspace_subdivision_matrix (hspace, hmsh);

% Fill the information for the boundaries
if (boundary)% && hmsh.ndim > 1)
  Nf = cumsum ([0, hspace.ndof_per_level]);
  for iside = 1:2*hmsh.ndim
    if (hmsh.ndim > 1)
      FtR_boundary = cell (size (FtR));
      for lev = 1:numel (FtR)
        [~,~,FtR_boundary{lev}] = intersect (FtR{lev}, hspace.space_of_level(lev).boundary(iside).dofs);
      end
      cells_boundary = cell (size (removed_cells));
      for lev = 1:numel (removed_cells)
        cells_boundary{lev} = get_boundary_indices (iside, hmsh.mesh_of_level(lev).nel_dir, removed_cells{lev});
      end
      hspace.boundary(iside) = hspace_coarsen2 (hspace.boundary(iside), hmsh.boundary(iside), FtR_boundary, cells_boundary);

      nlevels_aux = hspace.boundary(iside).nlevels;
    elseif (hmsh.ndim == 1)
      nlevels_aux = hspace.nlevels;
    end
      
    dofs = [];
    for lev = 1:nlevels_aux;
      [~,iact] = intersect (hspace.active{lev}, hspace.space_of_level(lev).boundary(iside).dofs);
      dofs = union (dofs, Nf(lev) + iact);
    end
    hspace.boundary(iside).dofs = dofs;
  end
  
else
  hspace.boundary = [];
end

% if (hspace.nlevels > hmsh.nlevels)
%   hspace = hspace_remove_empty_level (hspace, hmsh);
% end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function hspace = update_active_functions (hspace, marked_fun)
%
% This function updates the active dofs (hspace.active), their coefficients (hspace.coeff_pou) and deactivated dofs (hspace.deactivated) in each level when
% removing the functions in marked_fun. This function also updates hspace.nlevels, hspace.ndof and hspace.ndof_per_level
%
% Input:    hspace:    the fine space, an object of the class hierarchical_space
%           hmsh:      an object of the class hierarchical_mesh, already coarsened
%           marked_fun{lev}: indices of active functions of level lev to be removed
%
% Output:   hspace:    the coarsened space, an object of the class hierarchical_space
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
  if (~isempty (removed_cells{lev}))
    removed_funs = sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), removed_cells{lev});
    active{lev} = setdiff (active{lev}, removed_funs);
    active{lev-1} = union (active{lev-1}, funs_to_reactivate{lev-1});
    deactivated{lev-1} = setdiff (deactivated{lev-1}, funs_to_reactivate{lev-1});
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
