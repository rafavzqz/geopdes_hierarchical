% HIERARCHICAL_MESH: constructor of the class for hierarchical meshes.
%
%    hmsh = hierarchical_mesh (msh, geometry, [nsub=2*ones(1,msh.ndim)])
%
% INPUT:
%
%    msh:      the coarsest mesh (level 1), an object of the msh_structured class (see msh_structured)
%    geometry: an object of the geometry class (see geo_load)
%    nsub:     number of subdivisions between two different levels (by default 2)
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
%    active         (1 x nlevels cell-array) List of active elements on each level
%    deactivated    (1 x nlevels cell-array) List of removed cells on each level
%    msh_lev        (nlevels x 1 cell-array) msh_lev{ilev} is a structure
%    boundary       (2 x ndim array)    a hierarchical mesh representing the mesh on the boundary
%    geometry                           a copy of the geometry object (not so nice)
%
%    METHOD NAME
%    hmsh_plot_cells:       plot the hierarchical mesh (not efficient)
%    hmsh_refine:           refine the hierarchical mesh, given the set of cells to refine
%    hmsh_get_element_size: compute the size of the elements
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

function hmsh = hierarchical_mesh (msh, geometry, nsub)

if (nargin < 3)
  nsub = 2 * ones (1, msh.ndim);
end

hmsh.ndim = msh.ndim;
hmsh.rdim = msh.rdim;
hmsh.nlevels = 1;
hmsh.nsub = nsub;
hmsh.mesh_of_level = msh;
hmsh.nel = msh.nel;
hmsh.nel_per_level = [msh.nel];

aux = cell (hmsh.ndim, 1);
[aux{:}] = ind2sub ([msh.nel_dir, 1], 1:msh.nel); % The extra 1 makes it work in any dimension
hmsh.active{1} = (1:msh.nel)';
hmsh.deactivated{1} = [];
hmsh.msh_lev{1} = msh_evaluate_element_list (hmsh.mesh_of_level(1), hmsh.active{1});

hmsh.geometry = geometry; % CAN WE AVOID THIS? XXXXXX

if (~isempty (msh.boundary))
  if (msh.ndim > 1)
    for iside = 1:2*msh.ndim
      %%    ind =[2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind = [2 2 1 1] in 2D;
      ind = setdiff (1:hmsh.ndim, ceil (iside/2));
      hmsh.boundary(iside) = hierarchical_mesh (msh.boundary(iside), geometry.boundary(iside), nsub(ind));
    end
  elseif (msh.ndim == 1)
    hmsh.boundary(1).ndim = 0;
    hmsh.boundary(2).ndim = 0;
  end
else
  hmsh.boundary = [];
end

% Remove redundant information
hmsh.mesh_of_level(1).boundary = [];

hmsh = class (hmsh, 'hierarchical_mesh');

end