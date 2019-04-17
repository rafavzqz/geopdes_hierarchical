function [new_marked] = hrefine_rec(hmsh, hspace, Q_ind, lev_Q, m)

% HREFINE_REC: recursive refinement algorithm to guarantee admissible meshes.
%
%   [new_marked] = hrefine_rec (hmsh, hspace, Q_ind, lev_Q, m)
%
% INPUT:
%
%   hmsh:   object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the coarse space of hierarchical splines (see hierarchical_space)
%   Q_ind:  indices of the elements to be refined, all of the same level
%   Q_lev:  level of the elements in Q_ind
%   m:      admissibility class of the refined mesh
%
% OUTPUT:
%
%   new_marked: cell array with the indices, of the marked elements from coarser levels
%                to guarantee admissibility of class m with respect to Q_ind
%           
%  Algorithm from the paper: 
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

new_marked = cell (1,hmsh.nlevels);
lev_i=lev_Q-m+1;
neighbors = get_neighborhood (hmsh,hspace, Q_ind, lev_Q, m);
if numel(neighbors)>0
    new_new_marked = hrefine_rec (hmsh, hspace, neighbors, lev_i, m);
    for k=1:lev_i
        new_marked{k} = union (new_marked{k}, new_new_marked{k});
    end
end

if (numel(intersect(Q_ind,hmsh.active{lev_Q})) > 0)
    new_marked{lev_Q} = union (new_marked{lev_Q}, Q_ind);
end

end