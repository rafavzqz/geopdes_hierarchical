function [marked] = hrefine (hmsh, hspace, toberef, m)

% HREFINE: refinement algorithm to guarantee admissible meshes.
%
%   [marked_adm] = hrefine (hmsh, hspace, marked, m)
%
% INPUT:
%
%   hmsh:   object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the coarse space of hierarchical splines (see hierarchical_space)
%   marked: cell array with the indices, in the tensor product space, of the marked elements for each level
%   m:      admissibility class of the refined mesh
%
% OUTPUT:
%
%   marked_adm: cell array with the indices, in the tensor product space, of the marked elements for each level
%                such that the final mesh is admissible of class m
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

if m<2
    error('class of admissibility m<2')
end

marked = cell (1,hmsh.nlevels);
for l=1:hmsh.nlevels
    Q_ind = intersect (toberef{l},hmsh.active{l});
    if numel(Q_ind)>0
        new_marked = hrefine_rec (hmsh, hspace, Q_ind, l, m);
        for k=1:l
            marked{k} = union (marked{k}, new_marked{k});
        end
    end
end
end