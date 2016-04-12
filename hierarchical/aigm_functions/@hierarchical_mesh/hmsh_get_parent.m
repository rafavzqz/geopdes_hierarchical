% HMSH_GET_PARENT: compute the parent cell of a given cell of a certain level.
%
%     [parents, flag] = hmsh_get_parents (hmsh, lev, ind)
%
% Get the parent cell of a given cell of level lev, with the subdivision given by hmsh.nsub.
%  If several cells are given as input, return as output all their parent cells.
%
% INPUT:
%
%     hmsh: the hierarchical mesh (see hierarchical_mesh)
%     lev:  level of the subdivided cells to compute their parents
%     ind:  indices of the cells in the Cartesian grid
%
% OUTPUT:
%
%     parents: indices of the children, with the numbering of the Cartesian grid
%     flag:    a flag to tell whether all the input cells are active (1) 
%               active or deactivated (2), or if there is any passive cell (0)
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

function [parent, flag] = hmsh_get_parent (hmsh, lev, ind)

if (lev < 2 || lev > hmsh.nlevels)
  error ('The level should be between 2 and the number of levels of the mesh')
end

z = cell (hmsh.ndim, 1);
cells_sub = cell (hmsh.ndim, 1);
[cells_sub{:}] = ind2sub ([hmsh.mesh_of_level(lev).nel_dir, 1], ind); % The extra 1 makes it work in any dimension

parent = [];
for ii = 1:numel(cells_sub{1})
  aux = cell (hmsh.ndim, 1);
  for idim = 1:hmsh.ndim
      aux{idim} = floor ((cells_sub{idim} + hmsh.nsub(idim) - 1) / hmsh.nsub(idim));
  end
  [z{1:hmsh.ndim}] = ndgrid (aux{:});
  auxI = sub2ind ([hmsh.mesh_of_level(lev-1).nel_dir, 1], z{:});
  parent = union (parent, auxI(:));
end

if (nargout == 2)
  flag = all (ismember (ind, hmsh.active{lev}));
  if (~flag)
    flag = 2 *all (ismember (ind, union (hmsh.active{lev}, hmsh.deactivated{lev})));
  end
end

end
