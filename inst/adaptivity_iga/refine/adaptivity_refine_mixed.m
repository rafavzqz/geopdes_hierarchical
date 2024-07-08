% ADAPTIVITY_REFINE_MIXED: refine the hierarchical mesh and spaces, updating the corresponding structures hmsh and hspace.
%  Two or more spaces are refined at the same time.
%  The refinement can be done marking either elements or basis functions.
%  If refining by functions, the first space is considered.
%
%   [hmsh, hspaces] = adaptivity_refine (hmsh, hspaces, marked, adaptivity_data)
%
% INPUT:
%
%   hmsh:    object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   hspaces: cell-array with objects representing the coarse spaces of hierarchical splines (see hierarchical_space)
%   marked: cell array with the indices, in the tensor product space, of the marked elements/functions
%            for each level
%   adaptivity_data: a structure with the data for the adaptivity method.
%                    It contains the field 'flag', that can take the value
%                    'elements' or 'functions', depending on the refinement strategy.
%
% OUTPUT:
%
%   hmsh:    object representing the refined hierarchical mesh (see hierarchical_mesh)
%   hspaces: cell-array with objects representing the refined space of hierarchical splines (see hierarchical_space)
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2024 Rafael Vazquez
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

function [hmsh, hspaces] = adaptivity_refine_mixed (hmsh, hspaces, marked, adaptivity_data)

if (~iscell (hspaces))
  [hmsh, hspaces{1}] = adaptivity_refine (hmsh, hspaces, marked, adaptivity_data);
  return
elseif (numel(hspaces) == 1)
  [hmsh, hspaces{1}] = adaptivity_refine (hmsh, hspaces{1}, marked, adaptivity_data);
  return
end

switch (adaptivity_data.flag)
  case 'functions'
    marked_elements = compute_cells_to_refine (hspaces{1}, hmsh, marked);
  case 'elements'
    marked_elements = marked;
end

if (isfield (adaptivity_data, 'adm_class'))
  for ispace = 1:numel(hspaces)
    marked_elements = mark_admissible (hmsh, hspaces{ispace}, marked_elements, adaptivity_data);
  end
end
[hmsh, new_cells] = hmsh_refine (hmsh, marked_elements);

adaptivity_data.flag='elements';
for ispace = 1:numel(hspaces)
  marked_functions = compute_functions_to_deactivate (hmsh, hspaces{ispace}, marked_elements, 'elements');  
  hspaces{ispace} = hspace_refine (hspaces{ispace}, hmsh, marked_functions, new_cells);
end

end