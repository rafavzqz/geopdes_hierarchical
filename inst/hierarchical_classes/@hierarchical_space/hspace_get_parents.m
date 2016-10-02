% HSPACE_GET_PARENTS: compute the parents of a given set of functions of the same level.
%
%     [parents, flag] = hspace_get_parents (hspace, lev, ind)
%
% Get the parents of a given set of basis functions of level lev, with the
%  subdivision given by the "Proj" matrices.
% All the parent functions are stored in the same array.
%
% INPUT:
%
%     hspace: the hierarchical space (see hierarchical_space)
%     lev:    level of the children functions
%     ind:    indices of the functions in the tensor product space of level lev
%
% OUTPUT:
%
%     parents:  indices of the parents, with the numbering of the tensor product space
%     flag:     a flag to tell whether all the input functions are active (1) 
%               active or deactivated (2), or if there is any passive function (0)
%
% Copyright (C) 2016 Rafael Vazquez
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

function [parents, flag] = hspace_get_parents (hspace, lev, ind)

% % This computation would not work for the truncated basis
% Cmat = matrix_basis_change__ (hspace, lev);
% [~,parents] = find (Cmat(ind,:));
% parents = unique (parents);

ndim = size (hspace.Proj, 2);
z = cell (ndim, 1);
ind_sub = cell (ndim, 1);
[ind_sub{:}] = ind2sub ([hspace.space_of_level(lev).ndof_dir, 1], ind); % The extra 1 makes it work in any dimension

parents = [];
for ii = 1:numel(ind_sub{1})
  aux = cell (ndim, 1);
  for idim = 1:ndim
    aux{idim} = find (hspace.Proj{lev-1, idim}(ind_sub{idim}(ii),:));
  end
  [z{1:ndim}] = ndgrid (aux{:});
  auxI = sub2ind ([hspace.space_of_level(lev-1).ndof_dir, 1], z{:});
  parents = union (parents, auxI(:));
end

if (nargout == 2)
  flag = all (ismember (ind, hspace.active{lev}));
  if (~flag)
    flag = 2 *all (ismember (ind, union (hspace.active{lev}, hspace.deactivated{lev})));
  end
end

end
