% GET_CHILDREN_OF_CELLS: Given a list of elements, compute the position of their children in the 
%   full Cartesian grid of the next level.
%
%  = get_children_of_cells (hmsh, cell_numbers)
%
% INPUT:
%    hmsh:    XXXXXXXXXXXXXXXXXXXXXXXXXXXXxxXXXXXXXXXXXXXXXXXXXX
%    cell_numbers: numbering of the elements in the hierarchical mesh.
%
% OUTPUT:
%    indices: indices of the children, in the full Cartesian grid
%
% Copyright (C) 2015 Eduardo Garau, Rafael Vazquez
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



%%%%%%%
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% cell_numbers. Is this the globla numbering in the Cartesian grid?
% Or is it the position in hmsh.active{level}?

function get_children_of_cells (hmsh, level, cell_numbers)

nel_dir_coarse = hmsh.mesh_of_level(level).nel_dir;
nel_dir_fine = nel_dir_coarse .* hmsh.nsub;

% XXXXXXXX nsub could also depend on the level
indices = ind2sub ([nel_dir_coarse, 1], cell_numbers); % The extra 1 makes it work for any dimension
ind_fine = ((indices-1)*hmsh.nsub + 1):indices*hmsh.nsub;
sub2ind ([nel_dir_fine, 1], ind_fine);  % The extra 1 makes it work for any dimension

end
