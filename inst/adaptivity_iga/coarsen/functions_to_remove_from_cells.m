% FUNCTIONS_TO_REMOVE_FROM_CELLS: compute the indices of the active functions that have to be removed.
%
%   fun_indices = functions_to_remove_from_cells (hspace, hmsh, marked)
%
% INPUT:
%
%   hspace: object representing the (fine) space of hierarchical splines (see hierarchical_space)
%   hmsh:   object representing the coarsened hierarchical mesh (see hierarchical_mesh)
%   marked: cell array with the indices, in the Cartesian grid, of the elements that have to 
%             be removed for each level
%
% OUTPUT:
%   fun_indices: cell array with the indices of functions to be removed, for
%       each level, in the numbering of the tensor product space
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

function fun_indices = functions_to_remove_from_cells (hmsh, hspace, M)

fun_indices = cell (hspace.nlevels,1);

for lev = 2:hspace.nlevels
  if (~isempty(M{lev}))
    fun_indices{lev} = sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), M{lev});
    fun_indices{lev} = intersect (hspace.active{lev}, fun_indices{lev});
  end
end

end