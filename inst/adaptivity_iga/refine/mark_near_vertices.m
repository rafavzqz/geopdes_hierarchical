% MARK_NEAR_VERTICES: mark elements close to vertices to maintain the linear
%     independence of C^1 multipatch functions.
%
%   [marked, flag] = adaptivity_refine (hmsh, hspace, marked)
%
% INPUT:
%
%   hmsh:   object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the coarse space of hierarchical splines (see hierarchical_space)
%   marked: cell array with the indices, in the single level mesh, of the marked elements
%            for each level
%
% OUTPUT:
%
%   marked: cell array with the indices, in the single level mesh, of the marked elements
%            for each level
%   flag:   will return value 1 if any element had to be marked
%
% Copyright (C) 2021 Rafael Vazquez
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

function [marked, flag] = mark_near_vertices (hmsh, hspace, marked)

flag = false;
% Using brute force (check every vertex at every level)
vertices = hspace.space_of_level(1).vertices;
for ilev = 1:hmsh.nlevels
  new_marked = [];
  
  msh_lev = hmsh.mesh_of_level(ilev);
  [~, elems_adj_to_vertices] = msh_cells_near_vertex (msh_lev, vertices);
  
  marked_for_vertex = cellfun (@(x) intersect(marked{ilev}, x), elems_adj_to_vertices, 'UniformOutput', false);
  
  for ivert = 1:numel(vertices)
    if (~isempty (marked_for_vertex{ivert}))
      cells_on_patch = sp_get_vertex_neighbors (hspace.space_of_level(ilev), msh_lev, ivert);
      add_patch_elems = cellfun (@(x) any (ismember (marked_for_vertex{ivert}, x)), cells_on_patch);
      new_marked = union (new_marked, [cells_on_patch{add_patch_elems}]);
      flag = true;
    end
  end
  marked{ilev} = intersect (union (marked{ilev}, new_marked), hmsh.active{ilev});
end

end