% HSPACE_GET_CHILDREN: compute the children of a given set of functions of the same level.
%
%     [children, flag, children_of_function] = hspace_get_children (hspace, lev, ind)
%
% Get the children of a given set of basis functions of level lev, with the
%  subdivision given by the "Proj" matrices
% All the children functions are stored in the same array.
%
% INPUT:
%
%     hspace: the hierarchical space (see hierarchical_space_mp)
%     lev:    level of the functions to refine
%     ind:    indices of the functions in the multipatch space of level lev
%
% OUTPUT:
%
%     children: indices of the children, with the numbering of the tensor product space
%     flag:     a flag to tell whether all the input functions are active (1) 
%               active or deactivated (2), or if there is any passive function (0)
%     children_of_function: cell array with the children of each function
%
% Copyright (C) 2015, 2016, 2017 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2018 Cesare Bracco, Rafael Vazquez
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

function [children, flag, children_of_function] = hspace_get_children (hspace, lev, ind)


children = [];

children_of_function = cell (numel(ind), 1);

    ref_matrix = matrix_basis_change__ (hspace, lev+1);
    for ii = 1:numel(ind)
        auxI = find(ref_matrix(:,ind(ii)));
        children = union (children, auxI);
        children_of_function{ii} = auxI;  
    end
    


if (nargout == 2)
  flag = all (ismember (ind, hspace.active{lev}));
  if (~flag)
    flag = 2 *all (ismember (ind, union (hspace.active{lev}, hspace.deactivated{lev})));
  end
end

end
