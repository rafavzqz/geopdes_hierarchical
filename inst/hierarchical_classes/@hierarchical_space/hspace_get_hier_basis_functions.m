% HSPACE_GET_CHILDREN: compute the hierarchical active functions.
%
%     [hier_cells] = hspace_get_function_support (hspace, hmsh, active, lev, ind)
%
% Get the children of a given set of basis functions of level lev, with the
%  subdivision given by the "Proj" matrices
% All the children functions are stored in the same array.
%
% INPUT:
%
%     hspace: the hierarchical space (see hierarchical_space)
%     hmsh:   the hierarchical mesh (see hierarchical_mesh)
%     active: cell array of active function per level
%     lev:    level of the functions to refine
%     ind:    indices of the cell in the tensor product space of level lev
%
% OUTPUT:
%
%     hier_children: cell array with indices of the functions 
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

function [hier_basis] = hspace_get_hier_basis_functions (hspace, hmsh, active, lev, ind)

nactive_levels = numel(active) - numel(find(cellfun(@isempty, active)));
hier_basis = cell(1, nactive_levels);
basis_funct = sp_get_basis_functions(hspace.space_of_level(lev), hmsh.mesh_of_level(lev), ind);
hier_basis{lev} = intersect(basis_funct, active{lev});

for iLev = lev:-1:2
    support = hmsh_get_parent(hmsh, iLev, ind);
    basis_funct = sp_get_basis_functions(hspace.space_of_level(iLev-1), hmsh.mesh_of_level(iLev-1), support);
    
    hier_basis{iLev-1} = intersect(basis_funct, active{iLev-1});
    ind = support;
end

end
