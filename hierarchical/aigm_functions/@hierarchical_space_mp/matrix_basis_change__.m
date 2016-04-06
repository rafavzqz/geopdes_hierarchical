% MATRIX_BASIS_CHANGE__: compute the subdivision matrix between two consecutive levels.
%        This method is intended to remain private.
%
% function C = matrix_basis_change__ (hspace, hmsh, lev)
%
% Compute the new matrices to represent functions of level "lev"
% as linear combinations of splines (active and inactive) of the current level 
%
% INPUT:  
%
%   hspace: an object of the class hierarchical_space
%   hmsh:   an object of the class hierarchical_mesh
%   lev:    the level for which we compute the matrix
%
% OUTPUT:
%
%   C:    matrix to change basis from level lev-1 to level lev
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
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
function C = matrix_basis_change__ (hspace, hmsh, lev)

C = sparse (hspace.space_of_level(lev).ndof, hspace.space_of_level(lev-1).ndof);
for ipatch = 1:hmsh.npatch
  Cpatch = 1;
  Proj = hspace.Proj{lev-1, ipatch};
  for idim = 1:hmsh.ndim
    Cpatch = kron (Proj{idim}, Cpatch);
  end
  
  if (strcmpi (hspace.space_of_level(1).sp_patch{ipatch}.space_type, 'NURBS'))
    sp_patch = hspace.space_of_level(lev-1).sp_patch{ipatch};
    Wlev = spdiags (sp_patch.weights(:), 0, sp_patch.ndof, sp_patch.ndof);
    sp_patch = hspace.space_of_level(lev).sp_patch{ipatch};
    Wlev_fine = spdiags (sp_patch.weights(:), 0, sp_patch.ndof, sp_patch.ndof);
    Cpatch = Wlev_fine * Cpatch * Wlev;
  end
  C(hspace.space_of_level(lev).gnum{ipatch}, hspace.space_of_level(lev-1).gnum{ipatch}) = Cpatch;
end


if (hspace.truncated)
  indices = union (hspace.active{lev}, hspace.deactivated{lev});
  C(indices,:) = 0;
end

end
