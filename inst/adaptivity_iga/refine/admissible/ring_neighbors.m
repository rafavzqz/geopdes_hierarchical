%RING_NEIGHBORS evaluate a list of elements which have to be further
%refined to obtained a balanced function mesh with function support
%balancing parameter n_rb.
%
% the function implements Algorithm 7.2 of:
% Lorenzo et al., "Hierarchically refined and coarsened splines for moving
% interface problems, with particular application to phase field models of
% prostate tumor growth", 2017, CMAME pp. 1-46
%
%   balance_element_list = ring_neighbors( hmsh, ind, lev, n_rb )
%
%
% Copyright (C) 2017 Massimo Carraturo
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
function balance_element_list  = ring_neighbors( hmsh, ind, lev, n_rb )

% get parent's neighborhood
parent =  hmsh_get_parent (hmsh, lev, ind);
mhs = hmsh.mesh_of_level(lev-1);
res = [];
for k=1:numel(mhs.nel_dir)
res = [res n_rb];
end
[neigh, ~, ~] = neighbourND( parent, mhs.nel_dir, res );

balance_element_list =  hmsh_get_parent (hmsh, lev-1, neigh);

end

