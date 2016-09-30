% HMSH_REMOVE_EMPTY_LEVEL: remove the last level of the hierarchical mesh, only if it is empty.
%
%   hmsh = hmsh_remove_empty_level (hmsh)
%
% INPUT:
%
%   hmsh:    object representing the hierarchical mesh (see hierarchical_mesh)
%
% OUTPUT:
%
%   hmsh:   the object of the hierarchical mesh with one level less.
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

function hmsh = hmsh_remove_empty_level (hmsh)

% Remove last level if empty (this is not necessary)
if (isempty (hmsh.active{hmsh.nlevels}))
  for ibnd = 1:numel(hmsh.boundary)
    if (hmsh.boundary(ibnd).nel_per_level(hmsh.boundary(ibnd).nlevels) == 0)
       hmsh.boundary(ibnd) = hmsh_remove_empty_level (hmsh.boundary(ibnd));
    end
  end
    
  hmsh.nlevels = hmsh.nlevels - 1;
  hmsh.active = hmsh.active(1:hmsh.nlevels);
  hmsh.deactivated = hmsh.deactivated(1:hmsh.nlevels);
  hmsh.nel_per_level = hmsh.nel_per_level(1:hmsh.nlevels);
  hmsh.mesh_of_level = hmsh.mesh_of_level(1:hmsh.nlevels);
  hmsh.msh_lev = hmsh.msh_lev(1:hmsh.nlevels);  
end

end
