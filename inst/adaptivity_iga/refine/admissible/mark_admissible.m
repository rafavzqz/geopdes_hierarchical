% MARK_ADMISSIBLE: marking algorithm to guarantee admissible meshes.
%
%   [marked_adm] = mark_admissible (hmsh, hspace, marked, adaptivity_data)
%
% INPUT:
%
%   hmsh:      object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   hspace:    object representing the coarse space of hierarchical splines (see hierarchical_space)
%   marked:    cell array with the indices, in the tensor product space, of the marked elements for each level
%   adaptivity_data: struct that contains the following two fields
%    adm_class: admissibility class of the refined mesh. If m<2, there is no additional refinement.
%    adm_type:  admissibility type, either 'T-admissible' (default) or 'H-admissible'.
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
% The algorithm is the same as Algorithms 5-6 in the paper
%      C. Bracco, C. Giannelli, R. Vazquez
%      Refinement algorithms for adaptive isogeometric methods with
%      hierarchical splines, Axioms, 2018
% with a loop from the fine to the coarse level, to avoid recursive calls
% to the function.
%
% Copyright (C) 2017, 2018, 2019 Cesare Bracco, Rafael Vazquez
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

function [marked] = mark_admissible (hmsh, hspace, marked, adaptivity_data)

if (~isfield(adaptivity_data,'adm_class') || adaptivity_data.adm_class < 2)
  return
else
  m = adaptivity_data.adm_class;  
end

if (~isfield(adaptivity_data,'adm_type') || isempty (adaptivity_data.adm_type))
  adm_type = 'T-admissible';
else
  adm_type = adaptivity_data.adm_type;
end

if (~isa(hspace,'hierarchical_space_mp_C1'))
% Standard case
  for lev = hmsh.nlevels:-1:1
    neighbors = get_neighborhood (hmsh, hspace, marked{lev}, lev, m, adm_type);
    if (numel(neighbors) > 0)
      lev_m = lev - m + 1;
      new_marked = intersect (neighbors, hmsh.active{lev_m});
      marked{lev_m} = union (marked{lev_m}, new_marked);
    end
  end
else
% C1 multipatch. Start marking elements adjacent to a vertex, and then neighborhood for admissibility.
  for lev = hmsh.nlevels:-1:1
    new_marked=[];
    msh_lev = hmsh.mesh_of_level(lev);
    vertices = hspace.space_of_level(lev).vertices;
    [~, elems_adj_to_vertices] = msh_cells_near_vertex (msh_lev, vertices);

    marked_for_vertex = cellfun (@(x) intersect(marked{lev}, x), elems_adj_to_vertices, 'UniformOutput', false);

    for ivert = 1:numel(vertices)
      if (~isempty (marked_for_vertex{ivert}))
        cells_on_patch = sp_get_vertex_neighbors (hspace.space_of_level(lev), msh_lev, ivert);
        add_patch_elems = cellfun (@(x) any (ismember (marked_for_vertex{ivert}, x)), cells_on_patch);
        new_marked = union (new_marked, [cells_on_patch{add_patch_elems}]);
      end
    end
    marked{lev} = union (marked{lev}, intersect (new_marked, hmsh.active{lev}));
    
    neighbors = get_neighborhood (hmsh, hspace, marked{lev}, lev, m, adm_type);
    if (numel(neighbors) > 0)
      lev_m = lev - m + 1;
      new_marked = intersect (neighbors, hmsh.active{lev_m});
      marked{lev_m} = union (marked{lev_m}, new_marked);
    end
  end
end

end
