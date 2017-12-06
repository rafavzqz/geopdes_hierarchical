% HSPACE_GET_CHILDREN: compute the hierarchical active cell supports of function.
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
%     ind:    indices of the functions in the tensor product space of level lev
%
% OUTPUT:
%
%     hier_children: cell array with indices of the children per each
%     element
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

function [hier_cells] = hspace_get_function_support (hspace, hmsh, active, lev, ind)

nactive_levels = numel(active);% - numel(find(cellfun(@isempty, active)));
hier_cells = cell(1, nactive_levels);
function_support = sp_get_cells (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), ind);

for iLev = lev:nactive_levels
    hier_cells{iLev-lev+1} = intersect(function_support, hmsh.active{iLev});
    
    if ~(length(find(ismember(function_support,hmsh.active{iLev})))==numel(function_support))
        function_support = hmsh_get_children(hmsh, iLev, function_support(~ismember(function_support,hmsh.active{iLev})));
    else
        break;
    end
    
end

end
