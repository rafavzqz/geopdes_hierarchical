% MATRIX_BASIS_CHANGE__: compute the subdivision matrix between two consecutive levels.
%        This method is intended to remain private.
%
% function C = matrix_basis_change__ (hspace, lev)
%
% Compute the new matrices to represent functions of level "lev-1"
% as linear combinations of splines (active and inactive) of level "lev"
%
% INPUT:  
%
%   hspace: an object of the class hierarchical_space
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

function C = matrix_basis_change__ (hspace, lev, ind_coarse)

ndim = size (hspace.Proj, 2);

if (nargin < 3)
  C = 1;
  for idim = 1:ndim
    C = kron (hspace.Proj{lev-1,idim}, C);
  end

elseif (nargin == 3)
  sub_coarse = cell (ndim, 1);
  [sub_coarse{:}] = ind2sub ([hspace.space_of_level(lev-1).ndof_dir, 1], ind_coarse);
  
  rows = zeros (prod (hspace.space_of_level(lev).degree+1)*numel(ind_coarse), 1); cols = rows; vals = rows;
  ncounter = 0;
  for ii = 1:numel(ind_coarse)
    Caux = 1;
    for idim = 1:ndim
      Caux = kron (hspace.Proj{lev-1,idim}(:,sub_coarse{idim}(ii)), Caux);
    end
    [ir, ic, iv] = find (Caux);
    rows(ncounter+(1:numel(ir))) = ir;
    cols(ncounter+(1:numel(ir))) = ind_coarse(ii);
    vals(ncounter+(1:numel(ir))) = iv;
    ncounter = ncounter + numel (ir);
  end
  rows = rows(1:ncounter); cols = cols(1:ncounter); vals = vals(1:ncounter);
  C = sparse (rows, cols, vals, hspace.space_of_level(lev).ndof, hspace.space_of_level(lev-1).ndof);
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
