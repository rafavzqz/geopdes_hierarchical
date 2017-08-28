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

function [deact_marked, num] = active2deactivated_marking (marked, hmsh, hspace, adaptivity_data)

deact_marked = cell (hmsh.nlevels, 1);

for lev = 1:hmsh.nlevels-1
  if (~isempty(marked{lev+1}))
    switch (lower (adaptivity_data.flag))
      case 'elements'
        [parents, flag] = hmsh_get_parent (hmsh, lev+1, marked{lev+1});
        if (flag ~= 1)
          error ('Some nonactive elements were marked.')
        end
        
        parents = intersect (parents, hmsh.deactivated{lev});
        for ii = 1:numel(parents)
          if (all (ismember (hmsh_get_children (hmsh, lev, parents(ii)), marked{lev+1})))
            deact_marked{lev} = union(deact_marked{lev}, parents(ii));
          end
        end

      case 'functions'
        [parents, flag] = hspace_get_parents (hspace, lev+1, marked{lev+1});
        if (flag ~= 1)
          error ('Some nonactive functions were marked.')
        end

        parents = intersect (parents, hspace.deactivated{lev});
        for ii = 1:numel(parents)
          if (any (ismember (hspace_get_children (hspace, lev, parents(ii)), marked{lev+1})))
            deact_marked{lev} = union(deact_marked{lev}, parents(ii));
          end
        end
    end
  end
    
end
num = numel(cell2mat(deact_marked'));