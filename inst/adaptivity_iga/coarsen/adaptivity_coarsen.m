% ADAPTIVITY_COARSEN: coarsen the hierarchical mesh and space, updating the corresponding structures hmsh and hspace.
%  The refinement can be done marking either elements or basis functions.
%
%   [hmsh, hspace] = adaptivity_coarsen (hmsh, hspace, marked, adaptivity_data)
%
% INPUT:
%
%   hmsh:   object representing the fine hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the fine space of hierarchical splines (see hierarchical_space)
%   marked: cell array with the indices, in the tensor product space, of the marked elements/functions
%            to reactivate for each level
%   adaptivity_data: a structure with the data for the adaptivity method.
%                    In particular, it contains the field 'flag', that can take the value
%                    'elements' or 'functions', depending on the coarsening strategy.
%
% OUTPUT:
%
%   hmsh:   object representing the coarsened hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the coarsened space of hierarchical splines (see hierarchical_space)
%
% Copyright (C) 2016 Eduardo M. Garau, Rafael Vazquez
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

function [hmsh, hspace] = adaptivity_coarsen (hmsh, hspace, marked, adaptivity_data)

switch (adaptivity_data.flag)
  case 'functions'
    marked_elements = compute_cells_to_reactivate (hspace, hmsh, marked);
  case 'elements'
    marked_elements = marked;
end

[hmsh, removed_cells] = hmsh_coarsen (hmsh, marked_elements);

reactivated_fun = functions_to_reactivate_from_cells (hmsh, hspace, marked_elements);
hspace = hspace_coarsen (hspace, hmsh, reactivated_fun, removed_cells);

hmsh = hmsh_remove_empty_levels (hmsh);
hspace = hspace_remove_empty_levels (hspace, hmsh);

end