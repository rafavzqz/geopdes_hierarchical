% COMPUTE_FUNCTIONS_TO_DEACTIVATE: compute the indices of the active functions that have to be deactivated.
%
%   fun_indices = compute_functions_to_deactivate (hspace, hmsh, marked, flag)
%
% INPUT:
%
%   hspace: object representing the space of hierarchical splines (see hierarchical_space)
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   marked: cell array with the indices, in the tensor product space, of the marked 
%             element/functions for each level
%   flag:   the refinement strategy, marking either 'elements' or 'functions'
%
% OUTPUT:
%   fun_indices: cell array with the indices of functions to be deactivated, for
%       each level, in the numbering of the tensor product space
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

function I = compute_functions_to_deactivate (hmsh, hspace, M, flag)

I = cell (hspace.nlevels,1);

for i = 1:hspace.nlevels
  I{i} = zeros(0,1);
end

for lev = 1:hspace.nlevels
  if (~isempty(M{lev}))
  % Computation of candidate functions to be deactivated
    switch flag,
      case 'functions',
        I{lev} = sp_get_neighbors (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), M{lev});
        % Remove from I{lev} the nonactive functions and the active functions already selected for refinement
        I{lev} = setdiff (intersect (I{lev}, hspace.active{lev}), M{lev});
      case 'elements',
        I{lev} = sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), M{lev});
        % Remove from I{lev} the nonactive functions
        I{lev} = intersect (I{lev}, hspace.active{lev});
    end
    
    % Computation of functions that in fact have to be deactivated
    [dummy, cells_per_fun] = sp_get_cells (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), I{lev});
    flag_ell = cellfun (@(x) isempty (intersect (x, hmsh.active{lev})), cells_per_fun);
    I{lev} = I{lev} (flag_ell == 1);
        
    if strcmpi (flag,'functions')
      I{lev} = union (I{lev}, M{lev});
    end
  end
end

end