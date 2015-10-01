% HMSH_GET_ELEMENT_SIZE: compute the approximated size (length, area, volume) of each element.
%  The function uses the quadrature rule to compute the sizes.
%
%    [element_sizes, level_sizes] = hmsh_get_element_size (hmsh)
%
% INPUT:
%
%    hmsh:       object representing the hierarchical mesh (see hierarchical_mesh)
%
% OUTPUT:
%
%    element_sizes: the size of each element, for all the elements of the hierarchical mesh
%    level_sizes:   maximum size of cells of each level
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
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

function [element_sizes, level_sizes] = hmsh_get_element_size (hmsh)

element_sizes = zeros (hmsh.nel, 1);
level_sizes = zeros (hmsh.nlevels, 1);

Ne = cumsum ([0; hmsh.nel_per_level(:)]);
for lev = 1:hmsh.nlevels
  if (hmsh.nel_per_level(lev) > 0)
    ind_e = (Ne(lev)+1):Ne(lev+1);
    element_sizes(ind_e) = sum (hmsh.msh_lev{lev}.quad_weights .* hmsh.msh_lev{lev}.jacdet, 1);
    level_sizes(lev) = max (element_sizes(ind_e));
  end
end
