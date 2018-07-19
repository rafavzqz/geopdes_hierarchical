% ADAPTIVITY_REFINE_FSB: refine the hierarchical mesh and space, updating the corresponding structures hmsh and hspace.
%  The refinement can be done marking either elements or basis functions.
%
%   [hmsh, hspace, Cref] = adaptivity_refine_fsb (hmsh, hspace, marked, adaptivity_data)
%
% INPUT:
%
%   hmsh:   object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the coarse space of hierarchical splines (see hierarchical_space)
%   marked: cell array with the indices, in the tensor product space, of the marked elements/functions
%            for each level
%   adaptivity_data: a structure with the data for the adaptivity method.
%                    In particular, it contains the field 'flag', that can take the value
%                    'elements' or 'functions', depending on the refinement strategy.
%
% OUTPUT:
%
%   hmsh:   object representing the refined hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the refined space of hierarchical splines (see hierarchical_space)
%   Cref:   refinement matrix, to pass from the original space to the refined one
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
function [hmsh, hspace, Cref] = adaptivity_refine_fsb (hmsh, hspace, marked, adaptivity_data)

switch (adaptivity_data.flag)
    case 'functions'
        marked_elements = compute_cells_to_refine (hspace, hmsh, marked);
    case 'elements'
        marked_elements = marked;
end

[hmsh, new_cells] = hmsh_refine (hmsh, marked_elements);
balanced_marked = function_balancing (hmsh, hspace, adaptivity_data.adm);
[hmsh, new_cells_balanced] = hmsh_refine (hmsh, balanced_marked);
adm_mark = cell(numel(marked_elements), 1);
for iLevel = 1:numel(marked_elements)
    if (~isempty(balanced_marked{iLevel}) || ~isempty(marked_elements{iLevel}) )
        adm_mark{iLevel} = [marked_elements{iLevel}; balanced_marked{iLevel}];
    end
    if (~isempty(new_cells_balanced{iLevel}) || ~isempty(new_cells{iLevel}) )
        new_cells{iLevel} = [new_cells{iLevel}; new_cells_balanced{iLevel}];
    end

end

adaptivity_data.flag='elements';
marked_functions = compute_functions_to_deactivate (hmsh, hspace, adm_mark, adaptivity_data.flag);

if (nargout == 3)
    [hspace, Cref] = hspace_refine (hspace, hmsh, marked_functions, new_cells);
else
    hspace = hspace_refine (hspace, hmsh, marked_functions, new_cells);
end
end