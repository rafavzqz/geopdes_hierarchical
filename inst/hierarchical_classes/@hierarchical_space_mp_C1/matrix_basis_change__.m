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
%   hspace: an object of the class hierarchical_space_mp
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

% We need to change this function, to compute the coefficients in Mario's notes (refinement masks)
% It will write C^1 basis function of level lev-1, as linear combinations of C^1 basis functions of level lev.
% The output matrix must have size: ndof_{l} x ndof_{l-1}
% The coefficients are computed using the projection matrices for
%  univariate spaces in hspace.Proj, hspace.Proj0, hspace.Proj1
%    Proj  (hmsh.nlevels-1 x npatch cell-array) 
%           the coefficients relating 1D splines of two consecutive levels for each patch
%           Proj{l,i} is a cell-array of dimension ndim, with the information for
%           the univariate Projectors on the patch (see also hierarchical_space)
%    Proj0  Like Proj, for univariate splines of degree p, regularity r+1
%    Proj1  Like Proj, for univariate splines of degree p-1, regularity r


if (isa (hspace.space_of_level(1).sp_patch{1}, 'sp_scalar'))
  is_scalar = true;
  ndim = size (hspace.Proj{1}, 2);
% elseif (isa (hspace.space_of_level(1).sp_patch{1}, 'sp_vector'))
%   is_scalar = false;
%   ndim = size (hspace.Proj{1}, 2);
else
  error ('Unknown space type')
end

npatch = hspace.space_of_level(1).npatch;
ndim = size (hspace.Proj{1}, 2);

% For now we only use B-splines, and scalar spaces

% We have to change all this part following Mario's notes, and using Proj, Proj0 and Proj1
if (nargin < 3)
  C = sparse (hspace.space_of_level(lev).ndof, hspace.space_of_level(lev-1).ndof);
  for ipatch = 1:npatch
    Cpatch = 1;
    Proj = hspace.Proj{lev-1, ipatch};
    for idim = 1:ndim
      Cpatch = kron (Proj{idim}, Cpatch);
    end

    C(hspace.space_of_level(lev).gnum{ipatch}, hspace.space_of_level(lev-1).gnum{ipatch}) = Cpatch;
  end
  
  
% We have to change all this part following Mario's notes, and using Proj, Proj0 and Proj1
elseif (nargin == 3)
  nrows = hspace.space_of_level(lev).ndof; ncols = hspace.space_of_level(lev-1).ndof;
  size_alloc = numel (ind_coarse) * prod (hspace.space_of_level(lev).sp_patch{1}.degree + 1);
  rows = zeros (size_alloc, 1); cols = rows; vals = rows;
  ncounter = 0;
  for ipatch = 1:npatch
    nc_in_patch = ncounter+1;
    [ind_ptc, ~, local_indices] = intersect (ind_coarse, hspace.space_of_level(lev-1).gnum{ipatch});
    sub_coarse = cell (ndim, 1);
    [sub_coarse{:}] = ind2sub ([hspace.space_of_level(lev-1).sp_patch{ipatch}.ndof_dir, 1], local_indices);

    for ii = 1:numel(ind_ptc)
      Proj = hspace.Proj{lev-1, ipatch};
      Caux = 1;
      for idim = 1:ndim
        Caux = kron (Proj{idim}(:,sub_coarse{idim}(ii)), Caux);
      end
      
      [ir, ic, iv] = find (Caux);
%       rows(ncounter+(1:numel(ir))) = hspace.space_of_level(lev).gnum{ipatch}(ir);
%       cols(ncounter+(1:numel(ir))) = ic*ind_ptc(ii);
      rows(ncounter+(1:numel(ir))) = ir;
      cols(ncounter+(1:numel(ir))) = local_indices(ii);
      vals(ncounter+(1:numel(ir))) = iv;
      ncounter = ncounter + numel (ir);
    end

    rows(nc_in_patch:ncounter) = hspace.space_of_level(lev).gnum{ipatch}(rows(nc_in_patch:ncounter));
    cols(nc_in_patch:ncounter) = hspace.space_of_level(lev-1).gnum{ipatch}(cols(nc_in_patch:ncounter));
  end

  rows = rows(1:ncounter); cols = cols(1:ncounter); vals = vals(1:ncounter);
%   C = sparse (rows, cols, vals, hspace.space_of_level(lev).ndof, hspace.space_of_level(lev-1).ndof);
  if (~isempty (rows))
    C = accumarray ([rows,cols], vals.', [], @min, 0, true);
  else
    C = sparse (nrows, ncols);
  end
  if (size (C,1) < nrows || size(C,2) < ncols)
    C(nrows, ncols) = 0;
  end
end

if (hspace.truncated)
  indices = union (hspace.active{lev}, hspace.deactivated{lev});
  C(indices,:) = 0;
end

end
