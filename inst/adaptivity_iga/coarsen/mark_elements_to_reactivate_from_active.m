% ACTIVE2DEACTIVATED_MARKING: compute the deactivated entities (cells or functions) 
%  to be possibly reactivated, from a list of marked active entities.
%
%   [deact_marked, num] = active2deactivated_marking (marked, hmsh, hspace, adaptivity_data)
%
% INPUT:
%
%   marked:  cell-array with the indices of marked cells (or functions) for each level, in the tensor-product setting
%   hmsh:   object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the coarse space of hierarchical splines (see hierarchical_space)
%   adaptivity_data: a structure with the data for the adaptivity method.
%                    In particular, it contains the field:
%                     -'flag': elements or functions, according to marked
%
% OUTPUT:
%
%    deact_marked:  cell-array with the indices of marked cells (or functions) for each level to be reactivated, in the tensor-product setting
%    num         :  number of cells (or elements) to be reactivated
%
% This function computes, acording to adaptivity_data.flag, the deactivated entities to be possibly reactivated. 
% If flag = elements, all their children are (active and) marked. 
% If flag = functions, they have at least one child (active and) marked.  
%
% Copyright (C) 2016, 2017 Eduardo M. Garau, Rafael Vazquez
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
% admissible_to_reactivate = cell (hmsh.nlevels-1, 1);
% children_per_cell = cell (hmsh.nlevels-1, 1);

% % Elements that can be reactivated
% for lev = 1:hmsh.nlevels-1
%   [~,~,children_per_cell{lev}] = hmsh_get_children (hmsh, lev, hmsh.deactivated{lev});
%   ind = ~any (ismember (children_per_cell{lev}, hmsh.deactivated{lev+1}));
% %   ind = all (ismember (children_per_cell{lev}, hmsh.active{lev+1}));
%   admissible_to_reactivate{lev} = hmsh.deactivated{lev}(ind);
% end


% if (strcmpi (adaptivity_data.flag, 'elements'))
  for lev = 1:hmsh.nlevels-1
    if (~isempty(marked{lev+1}))
      [parents, flag] = hmsh_get_parent (hmsh, lev+1, marked{lev+1});
      if (flag ~= 1)
        error ('Some nonactive elements were marked.')
      end

      [~,~,children_per_cell] = hmsh_get_children (hmsh, lev, parents);    

      ind = all (ismember (children_per_cell, hmsh.active{lev+1}));
      parents = parents(ind);

      if (strcmpi (adaptivity_data.coarsening_flag, 'any'))
% To reactivate one cell, only one child needs to be marked
        deact_marked{lev} = parents;
      elseif (strcmpi (adaptivity_data.coarsening_flag, 'all'))
% To reactivate one cell, all the children must be marked
        ind2 = all (ismember (children_per_cell(:,ind), marked{lev+1}));
        deact_marked{lev} = parents(ind2);
      else
        error ('Unknown option for coarsening, in adaptivity_data.coarsening_flag')
      end
    end
  end
  
% elseif (strcmpi (adaptivity_data.flag, 'functions'))
%   for lev = 1:hmsh.nlevels-1
%     if (~isempty(marked{lev+1}))
%       [parents, flag] = hspace_get_parents (hspace, lev+1, marked{lev+1});
%       if (flag ~= 1)
%         error ('Some nonactive functions were marked.')
%       end
% 
%       parents = intersect (parents, hspace.deactivated{lev});
%       [~,~,children_per_function] = hspace_get_children (hspace, lev, parents);
%       
%       if (strcmpi (adaptivity_data.coarsening_flag, 'any'))
% % To reactivate one function, only one child needs to be marked
%         deact_marked{lev} = parents;
%       elseif (strcmpi (adaptivity_data.coarsening_flag, 'all'))
% % To reactivate one function, all the children must be marked
%         ind = cellfun (@(x) all (ismember (x, marked{lev+1})), children_per_function);
%         deact_marked{lev} = parents(ind);
%       else
%         error ('Unknown option for coarsening, in adaptivity_data.coarsening_flag')
%       end
%     end
%   end
% else
%   error ('Unknown option for coarsening, in adaptivity_data.flag')
% end

num = sum (cellfun (@numel, deact_marked));

% % deact_marked = cell (hmsh.nlevels, 1);
% 
% for lev = 1:hmsh.nlevels-1
%   [parents, flag] = hspace_get_parents (hspace, lev+1, marked{lev+1});
%   if (flag ~= 1)
%     error ('Some nonactive functions were marked.')
%   end
% 
%   parents = intersect (parents, hspace.deactivated{lev});
%   for ii = 1:numel(parents)
%     if (any (ismember (hspace_get_children (hspace, lev, parents(ii)), marked{lev+1})))
%       deact_marked{lev} = union(deact_marked{lev}, parents(ii));
%     end
%   end
    
end
