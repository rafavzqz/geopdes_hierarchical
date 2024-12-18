% HMSH_GET_CHILDREN: compute the children cells of a given set of cells of the same level.
%
%     [children, flag, children_of_cell] = hmsh_get_children (hmsh, lev, ind)
%
% Get the children cells of a given set of cells of level lev, with the
%  subdivision given by hmsh.nsub
%
% INPUT:
%
%     hmsh: the hierarchical mesh (see hierarchical_mesh_mp)
%     lev:  level of the cells to subdivide
%     ind:  indices of the cells in the global multipatch grid of that level
%
% OUTPUT:
%
%     children: indices of the children, with the numbering of the multipatch grid
%     flag:     a flag to tell whether all the input cells are active (1) 
%               active or deactivated (2), or if there is any passive cell (0)
%     children_of_cell: indices of the children (rows) for each cell in ind (column)
%
% Copyright (C) 2015, 2016, 2017 Eduardo M. Garau, Rafael Vazquez
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

function [children, flag, children_of_cell] = hmsh_get_children (hmsh, lev, ind)

if (lev < 1 || lev > hmsh.nlevels-1)
  error ('The level should be between 1 and the number of levels of the mesh minus one')
end

if (any (ind > hmsh.mesh_of_level(lev).nel))
  error ('There are some indices greater than the number of elements of the level')
end

children = [];
ndim = hmsh.ndim;
nsub = hmsh.nsub;

Nelem = cumsum ([0 hmsh.mesh_of_level(lev).nel_per_patch]);
Nelem_fine = cumsum ([0 hmsh.mesh_of_level(lev+1).nel_per_patch]);
z = cell (ndim, 1);
aux = cell (ndim, 1);
children_of_cell = zeros (prod (nsub), numel (ind));
for iptc = 1:hmsh.npatch
  [~,indices,position_on_ind] = intersect (Nelem(iptc)+1:Nelem(iptc+1), ind);

  cells_sub = cell (ndim, 1);
  [cells_sub{:}] = ind2sub ([hmsh.mesh_of_level(lev).msh_patch{iptc}.nel_dir, 1], indices); % The extra 1 makes it work in any dimension

  for ii = 1:numel(cells_sub{1})
    for idim = 1:ndim
      aux{idim} = nsub(idim)*(cells_sub{idim}(ii)-1)+1:nsub(idim)*(cells_sub{idim}(ii));
    end
    [z{1:ndim}] = ndgrid (aux{:});
    auxI = sub2ind ([hmsh.mesh_of_level(lev+1).msh_patch{iptc}.nel_dir, 1], z{:});
    children = union (children, auxI(:)+Nelem_fine(iptc));
    children_of_cell(:,position_on_ind(ii)) = auxI(:)' + Nelem_fine(iptc);
  end
end

if (nargout >= 2)
  flag = all (ismember (ind, hmsh.active{lev}));
  if (~flag)
    flag = 2 *all (ismember (ind, union (hmsh.active{lev}, hmsh.deactivated{lev})));
  end
end

end
