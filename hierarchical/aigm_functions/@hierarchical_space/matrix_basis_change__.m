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

function C = matrix_basis_change__ (hspace, lev)

ndim = size (hspace.Proj, 2);

C = 1;
for idim = 1:ndim
  C = kron (hspace.Proj{lev-1,idim}, C);
end

if (strcmpi (hspace.space_of_level(1).space_type, 'NURBS'))
  Wlev = spdiags (hspace.space_of_level(lev-1).weights(:), 0, hspace.space_of_level(lev-1).ndof, hspace.space_of_level(lev-1).ndof);
  Wlev_fine = spdiags (1./hspace.space_of_level(lev).weights(:), 0, hspace.space_of_level(lev).ndof, hspace.space_of_level(lev).ndof);
  C = Wlev_fine * C * Wlev;
end

if (hspace.truncated)
  indices = union (hspace.active{lev}, hspace.deactivated{lev});
  C(indices,:) = 0;
end

end