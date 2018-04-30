% MARK_ADMISSIBLE: marking algorithm to guarantee admissible meshes.
%
%   [marked_adm] = mark_admissible (hmsh, hspace, marked, m)
%
% INPUT:
%
%   hmsh:   object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the coarse space of hierarchical splines (see hierarchical_space)
%   marked: cell array with the indices, in the tensor product space, of the marked elements for each level
%   m:      admissibility class of the refined mesh. If m<2, there is no additional refinement.
%
% OUTPUT:
%
%   marked_adm: cell array with the indices, in the tensor product space, of the marked elements for each level
%                such that the final mesh is admissible of class m
%           
%  Adaption of the algorithm from the paper: 
%      A. Buffa and C. Giannelli
%      Adaptive isogeometric methods with hierarchical splines: error
%       estimator and convergence
%      Math. Models Meth. Appl. Sci., 2016
%
% Copyright (C) 2017, 2018 Cesare Bracco, Rafael Vazquez
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

function [marked] = mark_admissible (hmsh, hspace, marked, m)

if (m < 2)
  return
end

for lev = 1:hmsh.nlevels
  marked = mark_admissible_recursive (hmsh, hspace, marked, lev, m);
end

end


% MARK_ADMISSIBLE_recursive: recursive marking algorithm to guarantee admissible meshes.
%
%   [marked_adm] = mark_admissible_recursive (hmsh, hspace, marked, lev, m)
%
% INPUT:
%
%   hmsh:   object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the coarse space of hierarchical splines (see hierarchical_space)
%   marked: cell array with the indices, in the tensor product space, of the marked elements for each level
%   lev:    level of the cells for which we compute the neighborhood and support extension
%   m:      admissibility class of the refined mesh
%
% OUTPUT:
%
%   marked_adm: cell array with the indices, in the tensor product space, of the marked elements for each level
%                such that the final mesh is admissible of class m

function marked = mark_admissible_recursive (hmsh, hspace, marked, lev, m)
  neighbors = get_neighborhood (hmsh, hspace, marked{lev}, lev, m);
  if (numel(neighbors) > 0)
    lev_m = lev - m + 1;
    new_marked = intersect (neighbors, hmsh.active{lev_m});
    marked{lev_m} = union (marked{lev_m}, new_marked);
    marked = mark_admissible_recursive (hmsh, hspace, marked, lev_m, m);
  end
end