% HSPACE_REMOVE_EMPTY_LEVEL: remove the last level of the hierarchical space, only if it is empty.
%
%   hspace = hmsh_remove_empty_level (hspace, hmsh)
%
% INPUT:
%
%   hspace: object representing the hierarchical space (see hierarchical_space_mp)
%   hmsh:  object representing the hierarchical mesh, with the last level already removed (see hierarchical_mesh_mp)
%
% OUTPUT:
%
%   hspace: the object of the hierarchical space with one level less.
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

function hspace = hspace_remove_empty_level (hspace, hmsh)

if (hspace.nlevels == hmsh.nlevels+1 && isempty (hspace.active{end}))
  hspace.space_of_level(hspace.nlevels) = [];
  hspace.Proj(hspace.nlevels-1,:) = [];

  hspace.active(hspace.nlevels) = [];
  hspace.deactivated(hspace.nlevels) = [];
  hspace.ndof_per_level(hspace.nlevels) = [];
%   hspace.C(hspace.nlevels) = [];
  hspace.nlevels = hmsh.nlevels;
else
  warning ('The number of levels of the space in input should be one more than the number of levels of the mesh')
end

end
