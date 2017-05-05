% HSPACE_IN_FINER_MESH: construct a hierarchical space, in a mesh finer
%  than the hierarchical mesh that defines it.
%
%   hspace = hspace_in_finer_mesh (hspace, hmsh. hmsh_fine)
%
% INPUT:
%
%   hspace:    object representing the hierarchical space (see hierarchical_space)
%   hmsh:      object representing the hierarchical mesh on which the space is defined (see hierarchical_mesh)
%   hmsh_fine: object representing the finer hierarchical mesh where to evaluate the space (see hierarchical_mesh)
%
% OUTPUT:
%
%   hspace_on_fine: the object of the hierarchical space, with information
%    to be evaluated in the fine mesh
%
% The list of active and deactivated functions of the hierarchical space is unchanged.
% The difference is mainly in the matrix Csub, that is recomputed for the new mesh.
% Empty levels are added to the space, if needed.
% The function is intended to be used in particular situations, so the
%  boundary field is not computed, to save computational resources.
%
% Copyright (C) 2017 Rafael Vazquez
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

function hspace = hspace_in_finer_mesh (hspace, hmsh, hmsh_fine)

for ilev = 1:hmsh.nlevels
  if (~all (ismember (hmsh.active{ilev}, union (hmsh_fine.active{ilev}, hmsh_fine.deactivated{ilev}))))
    error ('The new mesh is not finer than the original one')
  end
end

% Adding empty levels
for ilev = hspace.nlevels+1:hmsh_fine.nlevels
  msh_level = hmsh_fine.mesh_of_level(ilev);
  degree = hspace.space_of_level(ilev-1).degree;
  [new_space, Proj] = sp_refine (hspace.space_of_level(ilev-1), msh_level, hmsh_fine.nsub, degree, degree-1);
  hspace.space_of_level(ilev) = new_space; clear new_space
  hspace.Proj(ilev-1,:) = Proj(:);
  
  hspace.nlevels = ilev;
  hspace.active{ilev} = [];
  hspace.deactivated{ilev} = [];
  hspace.ndof_per_level(ilev) = 0;
end

hspace.Csub = hspace_subdivision_matrix (hspace, hmsh_fine);

end