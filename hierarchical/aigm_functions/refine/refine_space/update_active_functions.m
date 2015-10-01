% function hspace = update_active_functions(hspace, hmsh, new_cells, marked_fun)
%
% This function updates the active dofs (hspace.active and hspace.globnum_active), their coefficients (hspace.coeff) and deactivated dofs (hspace.deactivated) in each level when
% refining the functions in marked_fun. This function also updates hspace.nlevels, hspace.ndof and hspace.ndof_per_level
%
% Input:    hspace:    the coarse mesh, an object of the class hierarchical_space
%           hmsh:      an object of the class hierarchical_mesh, already refined
%           new_cells: cells added to the refined mesh, see hmsh_refine
%           marked_fun{lev}: indices of active functions of level lev to be deactivated
%
% Output:   hspace:    the refined mesh, an object of the class hierarchical_space
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
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

function hspace = update_active_functions (hspace, hmsh, new_cells, marked_fun)

Nf = cumsum ([0; hspace.ndof_per_level(:)]);
W = cell (hspace.nlevels+1,1);

active = hspace.active;
deactivated = hspace.deactivated;

% El siguiente loop seguramente se puede evitar usando mat2cell
for lev = 1:hspace.nlevels
  ind_f = (Nf(lev)+1):Nf(lev+1);
  W{lev} = hspace.coeff_pou(ind_f);
  W{lev} = W{lev}(:);
end

for lev = 1:hspace.nlevels-1
  marked_fun{lev} = union (marked_fun{lev}, intersect (deactivated{lev}, active{lev}));
  if (~isempty(marked_fun{lev}))
    [dummy,indA] = ismember (marked_fun{lev}, active{lev});
%     if (~all (dummy))
%       disp('ERROR: update_active_functions: Some nonactive functions were selected');
%     end

    coefficients = 1;
    for idim = 1:hmsh.ndim
      coefficients = kron (hspace.Proj{lev,idim}, coefficients);
    end

    % Remove the corresponding functions from the active functions of level lev
    active{lev}(indA) = [];
    w = W{lev}(indA);
    W{lev}(indA) = [];
    deactivated{lev} = union (deactivated{lev}, marked_fun{lev});
    deactivated{lev} = deactivated{lev}(:);

    [ii,jj] = find (coefficients(:,marked_fun{lev}));
    children = arrayfun (@(x) ii(jj==x), 1:numel(marked_fun{lev}), 'UniformOutput', false);

% Add the children to the active functions, and compute the coefficients for the partition of unity
    for ifun = 1:numel (marked_fun{lev})
      Ichildren = children{ifun};
      c = coefficients(Ichildren, marked_fun{lev}(ifun));

      Ichildren_nonactive = setdiff (Ichildren, active{lev+1});
      if (~isempty (Ichildren_nonactive))
        nchildren_nonactive = numel (Ichildren_nonactive);
        active{lev+1} = vertcat (active{lev+1}, Ichildren_nonactive);
        W{lev+1} = vertcat (W{lev+1}, zeros(nchildren_nonactive,1));
      end
      [dummy, indices] = ismember (Ichildren, active{lev+1});
%       if (~all (dummy))
%         disp('ERROR: update_active_functions: Some nonactive functions were selected');
%       end
      W{lev+1}(indices) = W{lev+1}(indices) + w(ifun)*c;
    end % for ifun = 1:numel(marked_fun{lev})
  end %if
  
% For the classical hierarchical space, we activate functions of level lev+1, 
%  that are not children of any removed function of level lev
  if (strcmpi (hspace.type, 'standard') && ~isempty (new_cells{lev+1}))
        
    new_possible_active_fun = sp_get_basis_functions (hspace.space_of_level(lev+1), hmsh.mesh_of_level(lev+1), new_cells{lev+1});
%     new_possible_active_fun = setdiff (new_possible_active_fun, active{lev+1}, 'rows');
% XXXXXX I THINK THIS CHANGE (the union) WAS NECESSARY. ASK EDUARDO IF IT IS OK
    new_possible_active_fun = setdiff (new_possible_active_fun, union (active{lev+1}, deactivated{lev+1}));

    [dummy, elem] = sp_get_cells (hspace.space_of_level(lev+1), hmsh.mesh_of_level(lev+1), new_possible_active_fun);

    new_functions = cellfun (@(x) all (ismember (x, union (hmsh.active{lev+1}, hmsh.deactivated{lev+1}))), elem);
    active{lev+1} = vertcat (active{lev+1}, new_possible_active_fun(new_functions));
    W{lev+1} = vertcat (W{lev+1}, zeros (sum (new_functions), 1));
  end
  
  if (isempty (W{lev}))
    W{lev} = [];
    active{lev} = [];
  end
end % for lev


hspace.active = active(1:hspace.nlevels);
hspace.deactivated = deactivated(1:hspace.nlevels);
hspace.coeff_pou = cell2mat(W);
hspace.ndof_per_level = cellfun (@numel, hspace.active);
hspace.ndof = sum(hspace.ndof_per_level);

hspace.globnum_active = zeros (hspace.ndof,hmsh.ndim+1);
Nf = cumsum ([0; hspace.ndof_per_level(:)]);
for lev = 1:hspace.nlevels
  ind_f = (Nf(lev)+1):Nf(lev+1);
  if (~isempty(ind_f))
    hspace.globnum_active(ind_f, 1) = lev;
    globnum = cell (1,hmsh.ndim);
    [globnum{:}] = ind2sub (hspace.space_of_level(lev).ndof_dir, hspace.active{lev}(:));
    hspace.globnum_active(ind_f, 2:end) = cell2mat (globnum);
  end
end
