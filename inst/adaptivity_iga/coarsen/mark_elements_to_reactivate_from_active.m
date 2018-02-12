% MARK_ELEMENTS_TO_REACTIVATE_FROM_ACTIVE: compute the elements to
%     reactivate from the set of marked elements
%
%   [deact_marked, num] = mark_elements_to_reactivate_from_active (marked, hmsh, hspace, adaptivity_data)
%
% INPUT:
%
%   marked:  cell-array with the indices of marked cells (or functions) for each level, in the tensor-product setting
%   hmsh:   object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the coarse space of hierarchical splines (see hierarchical_space)
%   adaptivity_data: a structure with the data for the adaptivity method.
%
% OUTPUT:
%
%    deact_marked:  cell-array with the indices of marked cells (or functions) for each level to be reactivated, in the tensor-product setting
%    num         :  number of cells (or elements) to be reactivated
%
% Copyright (C) 2016, 2017, 2018 Eduardo M. Garau, Rafael Vazquez
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

function [deact_marked, num] = mark_elements_to_reactivate_from_active (marked, hmsh, hspace, adaptivity_data)

deact_marked = cell (hmsh.nlevels, 1);

for lev = hmsh.nlevels-1:-1:1
  if (~isempty(marked{lev+1}))
    [parents, flag] = hmsh_get_parent (hmsh, lev+1, marked{lev+1});
    if (flag ~= 1)
      error ('Some nonactive elements were marked.')
    end

    [~,~,children_per_cell] = hmsh_get_children (hmsh, lev, parents);    

    ind = all (ismember (children_per_cell, hmsh.active{lev+1}));
    parents = parents(ind);
    children_per_cell = children_per_cell(:,ind);

    if (strcmpi (adaptivity_data.coarsening_flag, 'any'))
% To reactivate one cell, only one child needs to be marked
      deact_marked{lev} = parents;
    elseif (strcmpi (adaptivity_data.coarsening_flag, 'all'))
% To reactivate one cell, all the children must be marked
      ind2 = all (ismember (children_per_cell, marked{lev+1}));
      parents = parents(ind2);
      children_per_cell = children_per_cell(:,ind2);
      deact_marked{lev} = parents;
    else
      error ('Unknown option for coarsening, in adaptivity_data.coarsening_flag')
    end
    
% Algorithm to recover admissible meshes
    if (isfield (adaptivity_data, 'adm') && adaptivity_data.adm > 1)
      lev_s = lev + adaptivity_data.adm;
      if (lev_s > hmsh.nlevels)
        continue
      else
        active_and_deact = union (hmsh.active{lev_s}, hmsh.deactivated{lev_s});
        supp_ext = support_extension (hmsh, hspace, children_per_cell(:), lev+1, lev+1);
        [~, descendants_of_cell] = hmsh_get_descendants (hmsh, supp_ext, lev+1, lev_s);
        keep_inds = [];
        for iel = 1:numel(deact_marked{lev})
          supp_ext_local = support_extension (hmsh, hspace, children_per_cell(:,iel), lev+1, lev+1);
          [~,ia,~] = intersect (supp_ext, supp_ext_local);
          if (isempty (intersect (descendants_of_cell(:,ia), active_and_deact)))
            keep_inds = [keep_inds, iel];
          end          
        end
        deact_marked{lev} = deact_marked{lev}(keep_inds);
      end
    end
    
  end
end
  
num = sum (cellfun (@numel, deact_marked));
    
end
