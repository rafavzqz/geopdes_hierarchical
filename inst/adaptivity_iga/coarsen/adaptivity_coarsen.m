% ADAPTIVITY_COARSEN: coarsen the hierarchical mesh and space, updating the corresponding structures hmsh and hspace.
%  The refinement can be done marking either elements or basis functions.
%
%   [hmsh, hspace, u] = adaptivity_coarsen (hmsh, hspace, marked, adaptivity_data)
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
%   u:      projected degrees of freedom
%
% Copyright (C) 2016 Eduardo M. Garau, Rafael Vazquez
%
% Copyright (C) 2017 Massimo Carraturo
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

function [hmsh, hspace, u_coarse] = adaptivity_coarsen (hmsh, hspace, marked, adaptivity_data)


switch (adaptivity_data.flag)
    case 'functions'
        [reactivated_fun, ~] = active2deactivated_marking(marked, hmsh, hspace, adaptivity_data);
        reactivated_elements = compute_cells_to_reactivate (hspace, hmsh, reactivated_fun);
    case 'elements'
        [reactivated_elements, ~] = active2deactivated_marking(marked, hmsh, hspace, adaptivity_data);
        reactivated_fun = functions_to_reactivate_from_cells (hmsh, hspace, reactivated_elements);
end

[hmsh, removed_cells] = hmsh_coarsen (hmsh, reactivated_elements);

if (nargout == 3)
    [hspace, u_coarse] = hspace_coarsen (hspace, hmsh, reactivated_fun, removed_cells, reactivated_elements);
else
    hspace = hspace_coarsen (hspace, hmsh, reactivated_fun, removed_cells, reactivated_elements);
end

hmsh = hmsh_remove_empty_levels (hmsh);
hspace = hspace_remove_empty_levels (hspace, hmsh);

end