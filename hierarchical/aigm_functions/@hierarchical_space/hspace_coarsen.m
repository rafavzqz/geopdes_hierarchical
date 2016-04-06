% HSPACE_COARSEN: coarsen the hierarchical space, updating the fields of the object.
%
%   hspace = hspace_coarsen (hspace, hmsh, functions_to_remove)
%
% INPUT:
%
%   hspace:         object representing the fine hierarchical space (see hierarchical_space)
%   hmsh:           object representing the hierarchical mesh, already coarsened (see hierarchical_mesh)
%   funs_to_remove: cell array with the indices, in the tensor product setting, of the functions to remove for each level
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
  Nf = cumsum ([0, hspace.ndof_per_level]);
  for iside = 1:2*hmsh.ndim
    if (hmsh.ndim > 1)
      M_boundary = cell (size (M));
      for lev = 1:numel (M)
        [~,~,M_boundary{lev}] = intersect (M{lev}, hspace.space_of_level(lev).boundary(iside).dofs);
      end
      hspace.boundary(iside) = hspace_coarsen (hspace.boundary(iside), hmsh.boundary(iside), M_boundary);

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

if (hspace.nlevels > hmsh.nlevels)
  hspace = hspace_remove_empty_level (hspace, hmsh);
end

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
% Input:  hspace: an object of the class hierarchical_space
%         hmsh:   an object of the class hierarchical_mesh
%         lev:    the level for which we compute the matrix
%
% Output:   C:    matrix to change basis from level lev-1 to level lev

function C = matrix_basis_change (hspace, hmsh, lev)

C = 1;
for idim = 1:hmsh.ndim
  C = kron (hspace.Proj{lev-1,idim}, C);
end

if (strcmpi (hspace.space_of_level(1).space_type, 'NURBS'))
  Wlev = spdiags (hspace.space_of_level(lev-1).weights(:), 0, hspace.space_of_level(lev-1).ndof, hspace.space_of_level(lev-1).ndof);
  Wlev_fine = spdiags (1./hspace.space_of_level(lev).weights(:), 0, hspace.space_of_level(lev).ndof, hspace.space_of_level(lev).ndof);
  C = Wlev_fine * C * Wlev;
end

if (hspace.truncated)
  indices = union (hspace.active{lev}, hspace.deactivated{lev});
  C(indices,:) = 0;
end

end