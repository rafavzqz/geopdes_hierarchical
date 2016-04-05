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

function hspace = hspace_coarsen (hspace, hmsh, M)

boundary = ~isempty (hspace.boundary);

% Update active functions
hspace = update_active_functions (hspace, hmsh, M);

% Update C, the matrices for changing basis
C = cell (hmsh.nlevels, 1);
C{1} = speye (hspace.space_of_level(1).ndof);
C{1} = C{1}(:,hspace.active{1});

for lev = 2:hmsh.nlevels
  I = speye (hspace.space_of_level(lev).ndof);
  aux = matrix_basis_change (hspace, hmsh, lev);
  C{lev} = [aux*C{lev-1}, I(:,hspace.active{lev})];
end
hspace.C = C;
clear C

% Fill the information for the boundaries
if (boundary)% && hmsh.ndim > 1)
  if (hmsh.ndim > 1)
    M_boundary = cell (size (M));
    levels = find (~cellfun (@isempty, M));
    for lev = levels(:).'
      [~,~,M_boundary{lev}] = intersect (M{lev}, hspace.space_of_level(lev).boundary.dofs);
    end
    
    hspace.boundary = hspace_coarsen (hspace.boundary, hmsh.boundary, M_boundary);

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

function hspace = update_active_functions (hspace, hmsh, marked_fun)

active = hspace.active;
deactivated = hspace.deactivated;

for lev = hspace.nlevels:-1:2
  active{lev} = setdiff (active{lev}, marked_fun{lev});
%   parents = hspace_get_parents (hspace, lev, marked_fun{lev});
  Cmat = matrix_basis_change (hspace, hmsh, lev);
  [~,parents] = find (Cmat(marked_fun{lev},:));
  parents = intersect (parents, deactivated{lev-1});
  active{lev-1} = union (active{lev-1}, parents);
  deactivated{lev-1} = setdiff (deactivated{lev-1}, parents);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function C = matrix_basis_change (hspace, hmsh, lev)
%
% Compute the new matrices to represent functions of the previous level 
% as linear combinations of splines (active and inactive) of the current level 
%
% Input:  hspace: an object of the class hierarchical_space_mp
%         hmsh:   an object of the class hierarchical_mesh_mp
%         lev:    the level for which we compute the matrix
%
% Output:   C:    matrix to change basis from level lev-1 to level lev

function C = matrix_basis_change (hspace, hmsh, lev)

C = sparse (hspace.space_of_level(lev).ndof, hspace.space_of_level(lev-1).ndof);
for ipatch = 1:hmsh.npatch
  Cpatch = 1;
  Proj = hspace.Proj{lev-1, ipatch};
  for idim = 1:hmsh.ndim
    Cpatch = kron (Proj{idim}, Cpatch);
  end
  
  if (strcmpi (hspace.space_of_level(1).sp_patch{ipatch}.space_type, 'NURBS'))
    sp_patch = hspace.space_of_level(lev-1).sp_patch{ipatch};
    Wlev = spdiags (sp_patch.weights(:), 0, sp_patch.ndof, sp_patch.ndof);
    sp_patch = hspace.space_of_level(lev).sp_patch{ipatch};
    Wlev_fine = spdiags (sp_patch.weights(:), 0, sp_patch.ndof, sp_patch.ndof);
    Cpatch = Wlev_fine * Cpatch * Wlev;
  end
  C(hspace.space_of_level(lev).gnum{ipatch}, hspace.space_of_level(lev-1).gnum{ipatch}) = Cpatch;
end


if (hspace.truncated)
  indices = union (hspace.active{lev}, hspace.deactivated{lev});
  C(indices,:) = 0;
end

end
