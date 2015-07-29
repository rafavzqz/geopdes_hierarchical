% HIERARCHICAL_MESH: constructor of the class for hierarchical meshes.
%
%    function hmsh = hierarchical_mesh (msh, geometry)
%
% INPUT:
%
%    msh:      the coarsest mesh (level 1), an object of the msh_structured class (see msh_structured)
%    geometry: an object of the geometry class (see geo_load)
%
% OUTPUT:
%  
%    hmsh: hierarchical_mesh object, which contains the following fields and methods
% 
%    FIELD_NAME     TYPE                DESCRIPTION
%    ndim           (scalar)            dimension of the parametric space
%    rdim           (scalar)            dimension of the physical space
%    nlevels        (scalar)            the number of levels
%    nsub           (1 x ndim array)    number of subdivisions between two different levels
%    mesh_of_level  (1 x nlevels mesh)  Cartesian mesh of each level (see msh_cartesian)
%    nel            (scalar)            total number of active cells  
%    nel_per_level  (1 x nlevels array) number of active cells on each level
%    globnum_active (nel x (dim+1))     global tensor-product numbering of active cells and their corresponding level
%    active         (1 x nlevels cell-array) List of active elements on each level
%    deactivated    (1 x nlevels cell-array) List of removed cells on each level
%    msh_lev        (nlevels x 1 cell-array) msh_lev{ilev} is a structure
%    geometry XXXXXXXXXXx
%    boundary       (2 x ndim array)    a hierarchical mesh representing the mesh on the boundary
%
%    METHOD NAME
%    hmsh_plot    XXXXXXXXX
%    hmsh_refine  XXXXXXXXX
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

function hmsh = hierarchical_mesh (msh, geometry)

hmsh.ndim = msh.ndim;
hmsh.rdim = msh.rdim;
hmsh.nlevels = 1;
hmsh.nsub = 2 * ones (1, hmsh.ndim); % Provisorio XXXXXX
hmsh.mesh_of_level = msh;
hmsh.nel = msh.nel;
hmsh.nel_per_level = [msh.nel];

aux = cell (hmsh.ndim, 1);
[aux{:}] = ind2sub ([msh.nel_dir, 1], 1:msh.nel); % The extra 1 makes it work in any dimension
hmsh.globnum_active = [ones(msh.nel, 1), cell2mat(aux)'];
hmsh.active{1} = (1:msh.nel)';
hmsh.deactivated{1} = [];
hmsh.msh_lev{1} = msh_evaluate_element_list (hmsh.mesh_of_level(1), hmsh.active{1});

hmsh.geometry = geometry; % CAN WE AVOID THIS? XXXXXX

if (~isempty (msh.boundary) && msh.ndim > 1)
  for iside = 1:2*msh.ndim
    hmsh.boundary(iside) = hierarchical_mesh (msh.boundary(iside), geometry.boundary(iside));
  end
else
  hmsh.boundary = [];
end

end