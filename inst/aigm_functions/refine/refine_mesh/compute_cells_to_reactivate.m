% COMPUTE_CELLS_TO_REACTIVATE: compute the indices of the cells that have to be reactivated, when 
%  the coarsening strategy is based on marking functions.
%
%   [marked_elem, indices] = compute_cells_to_reactivate (hspace, hmsh, marked)
%
% INPUT:
%
%   hspace: object representing the space of hierarchical splines (see hierarchical_space)
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   marked_fun: cell array with the indices, in the tensor product space, of the marked functions
%            to reactivate for each level
%
% OUTPUT:
%
%   marked_elem: cell array with the indices, in the Cartesian grid, of the elements
%                 to be reactivated for each level
%   indices:     relative position of the elements in the numbering of deactivated cells
%                 within the level, that is, in hmsh.deactivated{lev}
%
% Copyright (C) 2016 Eduardo M. Garau, Rafael Vazquez
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

function [ME, ind] = compute_cells_to_reactivate (hspace, hmsh, MF)

ME = cell (hmsh.nlevels, 1);
ind = cell (hmsh.nlevels, 1);

for lev = 1:hspace.nlevels-1
  if (~isempty(MF{lev}))
    elems = sp_get_cells (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), MF{lev});
  
    stay_deactivated = setdiff (hspace.deactivated{lev}, MF{lev});
    for iel = 1:numel(elems)
      children = hmsh_get_children (hmsh, lev, elems(iel));
      if (all (ismember (children, hmsh.active{lev+1}))) % if (isempty (intersect (children, hmsh.deactivated{lev+1}))) 
        funs = sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), elems(iel));
        if (~any (ismember (funs, stay_deactivated))) % if (isempty (intersect (funs, stay_deactivated))) 
          ME{lev} = vertcat (ME{lev}, elems(iel));
        end
      end
    end
    [~, ~, ind{lev}] = intersect (ME{lev}, hmsh.deactivated{lev});
  end
end

end