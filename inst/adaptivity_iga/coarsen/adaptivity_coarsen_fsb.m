%ADAPTIVITY_COARSEN_FSB coarse the hierarchical mesh and space, updating the corresponding structures hmsh and hspace.
%  The refinement can be done marking either elements or basis functions.
%
%   [hmsh, hspace, u] = adaptivity_coarsen_fsb (hmsh, hspace, marked, adaptivity_data)
%
% INPUT:
%
%   hmsh:   object representing the fine hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the fine space of hierarchical splines (see hierarchical_space)
%   marked: cell array with the indices, in the tensor product space, of the marked elements/functions
%            for each level
%   adaptivity_data: a structure with the data for the adaptivity method.
%                    In particular, it contains the field 'flag', that can take the value
%                    'elements' or 'functions', depending on the refinement strategy.
%
% OUTPUT:
%
%   hmsh:   object representing the coarsened hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the coarsened space of hierarchical splines (see hierarchical_space)
%   u_coarse:   coarsened dofs
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
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
function [hmsh, hspace, C_coarse] = adaptivity_coarsen_fsb (hmsh_fine, hspace, marked, adaptivity_data)

switch (adaptivity_data.flag)
    case 'functions'
        [reactivated_fun, ~] = active2deactivated_marking(marked, hmsh_fine, hspace, adaptivity_data);
        reactivated_elements = compute_cells_to_reactivate (hspace, hmsh_fine, reactivated_fun);
    case 'elements'
        [reactivated_elements, ~] = active2deactivated_marking(marked, hmsh_fine, hspace, adaptivity_data);
        reactivated_fun = functions_to_reactivate_from_cells (hmsh_fine, hspace, reactivated_elements);
end

[hmsh_coarse, removed_cells] = hmsh_coarsen (hmsh_fine, reactivated_elements);
switch(adaptivity_data.coarse_flag)
    case 'MS_all'
        [hspace_coarse] = hspace_coarsen_MS_all_active (hspace, hmsh_coarse, reactivated_fun, removed_cells);
    case 'MS_old'
        [hspace_coarse] = hspace_coarsen_MS_old_active (hspace, hmsh_coarse, reactivated_fun, removed_cells);
    case 'L2_global'
        hspace_coarse = hspace_coarsen (hspace, hmsh_coarse, reactivated_fun, removed_cells );
end

balanced_marked = function_balancing (hmsh_coarse, hspace_coarse, adaptivity_data.adm);
reactivated_balanced_elements = cell(numel(marked), 1);
for iLevel = 1:numel(marked)-1
    if ~isempty(reactivated_elements{iLevel})
         reactivated_balanced_elements{iLevel} = setdiff(reactivated_elements{iLevel}, balanced_marked{iLevel});
    end
end

adaptivity_data.flag='elements';
reactivated_balanced_fun = functions_to_reactivate_from_cells (hmsh_fine, hspace, reactivated_balanced_elements);

[hmsh, removed_cells] = hmsh_coarsen (hmsh_fine, reactivated_balanced_elements);

if (nargout == 3)
    switch(adaptivity_data.coarse_flag)
        case 'MS_all'
            [hspace, C_coarse] = hspace_coarsen_MS_all_active (hspace, hmsh, reactivated_balanced_fun, removed_cells);
        case 'MS_old'
            [hspace, C_coarse] = hspace_coarsen_MS_old_active (hspace, hmsh, reactivated_balanced_fun, removed_cells);
        case 'L2_global'
              hspace_fine = hspace;
              hspace = hspace_coarsen (hspace, hmsh, reactivated_balanced_fun, removed_cells );
              M = op_u_v_hier (hspace, hspace, hmsh);
              G = op_u_v_hier (hspace_fine, hspace_in_finer_mesh(hspace, hmsh, hmsh_fine), hmsh_fine);
              C_coarse = M \ G; C_coarse(abs(C_coarse) < 1e-12) = 0;
%               u_coarse =  Ccoar*hspace.dofs;
    end
else
    switch(adaptivity_data.coarse_flag)
        case 'MS_all'
            [hspace] = hspace_coarsen_MS_all_active (hspace, hmsh, reactivated_balanced_fun, removed_cells);
        case 'MS_old'
            [hspace] = hspace_coarsen_MS_old_active (hspace, hmsh, reactivated_balanced_fun, removed_cells);
        case 'L2_global'
            hspace = hspace_coarsen (hspace, hmsh, reactivated_balanced_fun, removed_cells );
    end
end

hmsh = hmsh_remove_empty_levels (hmsh);
hspace = hspace_remove_empty_levels (hspace, hmsh);

end

