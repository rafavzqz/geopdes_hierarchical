function [ list_of_cells ] = get_neighborhood (hmsh, hspace, Q_ind, lev_Q, m)

% GET_NEIGHBORHOOd: compute the neighborhood of the cells of a hierarchical mesh
%
%   [list_of_cells] = get_neighborhood (hmsh, hspace, Q_ind, Q_lev, m)
%
% INPUT:
%
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the space of hierarchical splines (see hierarchical_space)
%   Q_ind:  indices of the input elements, all of the same level
%   Q_lev:  level of the elements in Q_ind
%   m:      admissibility class of the hierarchical mesh
%
% OUTPUT:
%
%   list_of_cells: cell array with the indices, in the tensor product space,
%      of the active elements in the neighborhood of the elements in Q_ind
%           
%  The definition of the neighborhood is from the paper: 
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


lev_s = lev_Q-m+2;
if lev_s <= 1
    list_of_cells = [];
else
    all_el_lev = 1:hmsh.mesh_of_level(lev_s).nel;
    inactive_el = setdiff (all_el_lev, union (hmsh.active{lev_s}, hmsh.deactivated{lev_s}));
    el = intersect (inactive_el, support_extension(hmsh,hspace,Q_ind,lev_Q, lev_s));
    ancestors = hmsh_get_parent (hmsh, lev_s, el);
    list_of_cells = intersect (ancestors, hmsh.active{lev_s-1});
end

end