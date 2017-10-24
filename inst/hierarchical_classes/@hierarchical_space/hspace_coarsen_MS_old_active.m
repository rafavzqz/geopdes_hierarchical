% HSPACE_COARSEN: coarsen the hierarchical space, updating the fields of the object.
%
%   hspace = hspace_coarsen (hspace, hmsh, funs_to_reactivate, removed_cells)
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

% Here I implement the QI by Speleers and Manni
% The local projection is the L2 projection
% For functions that were active, I keep the same coefficient without changes
% For reactivated functions, I use as an interval the elements of the same
%  level within the support


function [hspace, Ccoarse] = hspace_coarsen_MS_old_active (hspace, hmsh, FtR, removed_cells)

boundary = ~isempty (hspace.boundary);

% Update active functions
if (nargout == 2)
  [hspace, Ccoarse] = update_active_functions (hspace, hmsh, FtR, removed_cells);
else
  hspace = update_active_functions (hspace, hmsh, FtR, removed_cells);
end

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
      hspace.boundary(iside) = hspace_coarsen_MS_old_active (hspace.boundary(iside), hmsh.boundary(iside), FtR_boundary, cells_boundary);

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

function [hspace, Ccoarse] = update_active_functions (hspace, hmsh, funs_to_reactivate, removed_cells)

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

ndof_old = sum (cellfun (@numel, hspace.active));
ndof = sum (cellfun (@numel, active));

if (nargout == 2)
  if (~hspace.truncated)
    error ('THB-splines are required to compute the coarsening matrix');
  end
    
  Ccoarse = sparse (ndof, ndof_old);
  for lev = hspace.nlevels:-1:1
    [~, inew, iold] = intersect (active{lev}, hspace.active{lev});
    ndof_until_lev_old = sum (hspace.ndof_per_level(1:lev-1));
    ndof_until_lev = sum (cellfun (@numel, active(1:lev-1)));
%     if (lev < hspace.nlevels)
%       ndof_until_next_lev = sum (hspace.ndof_per_level(1:lev+1));
%     end
% Functions already active keep the same coefficient
    Ccoarse(ndof_until_lev+inew,ndof_until_lev_old+iold) = speye (numel(inew), numel(inew));
    
% For other functions, I compute the integrals
   if (lev < hspace.nlevels)
    reactivated_cells = hmsh_get_parent (hmsh, lev+1, removed_cells{lev+1});
    old_active_cells = setdiff (hmsh.active{lev}, reactivated_cells);
  
    msh_lev1 = msh_evaluate_element_list (hmsh.mesh_of_level(lev+1), removed_cells{lev+1});
    
    ireact = setdiff (1:numel(active{lev}), inew);
    [all_cells, cells_per_function] = ...
      sp_get_cells (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), active{lev}(ireact));
    for ifun = 1:numel(ireact)
      cells = intersect (cells_per_function{ifun}, hmsh.active{lev});
      msh_support = msh_restrict_to_cells (hmsh.msh_lev{lev}, cells);
      sp_support = sp_evaluate_element_list (hspace.space_of_level(lev), msh_support);
      indices = unique (sp_support.connectivity);
      [~,position] = ismember (sp_support.connectivity, indices);
      sp_support.ndof = numel (indices);
      sp_support.connectivity = position;
      M_fun = op_u_v (sp_support, sp_support, msh_support, ones (msh_support.nqn, msh_support.nel));

% Computation in cells that have been reactivated. We pass to the
%  children, as they involve functions of level lev+1
      new_cells = setdiff (cells, old_active_cells);
      children = hmsh_get_children (hmsh, lev, new_cells);
      msh_children = msh_restrict_to_cells (msh_lev1, children);
      sp_lev1 = sp_evaluate_element_list (hspace.space_of_level(lev+1), msh_children);

% To use matrix_basis_change__, we must pass to the standard basis
%  otherwise some coefficients are set to zero
      hspace.truncated = false;
      Cproj = matrix_basis_change__ (hspace, lev+1, indices);
      Cproj = Cproj(:,indices);
      hspace.truncated = true;
    
      Gprev = op_u_v (sp_lev1, sp_lev1, msh_children, ones (msh_children.nqn, msh_children.nel));
      Gprev = Cproj.' * Gprev * hspace.Csub{lev+1};

      [~,jj,~] = find (Gprev);
      aux = M_fun \ Gprev(:,unique(jj));
      aux (abs(aux) < 1e-12) = 0;
      index = find (indices == active{lev}(ireact(ifun)));
      Ccoarse(ndof_until_lev+ireact(ifun), unique(jj)) = aux(index,:);
    end
   end
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
