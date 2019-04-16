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

function varargout = matrix_change_two_levels__ (sp_coarse, sp_fine, Proj, ind_coarse)

if (isa (sp_coarse, 'sp_scalar'))
  is_scalar = true;
elseif (isa (sp_coarse, 'sp_vector'))
  is_scalar = false;
  ncomp_param = sp_coarse.ncomp_param;
else
  error ('Unknown space type')
end
ndim = size (Proj, 2);


if (nargin < 4)
  if (is_scalar)
    C = 1;
    for idim = 1:ndim
      C = kron (Proj{1,idim}, C);
    end
    if (strcmpi (sp_coarse.space_type, 'NURBS'))
      Wlev = spdiags (sp_coarse.weights(:), 0, sp_coarse.ndof, sp_coarse.ndof);
      Wlev_fine = spdiags (1./sp_fine.weights(:), 0, sp_fine.ndof, sp_fine.ndof);
      C = Wlev_fine * C * Wlev;
    end
    
  else
    Caux = cell (ncomp_param, 1);
    for icomp = 1:ncomp_param
      spc_scalar = sp_coarse.scalar_spaces{icomp};
      spf_scalar = sp_fine.scalar_spaces{icomp};
      Caux{icomp} = matrix_change_two_levels__ (spc_scalar, spf_scalar, Proj(icomp,:));
    end
    C = blkdiag (Caux{:});
  end

  
  
elseif (nargin == 4)
  if (is_scalar)
    sub_coarse = cell (ndim, 1);
    [sub_coarse{:}] = ind2sub ([sp_coarse.ndof_dir, 1], ind_coarse);
  
    rows = zeros (prod (sp_fine.degree+1)*numel(ind_coarse), 1); cols = rows; vals = rows;
    ncounter = 0;
    for ii = 1:numel(ind_coarse)
      Caux = 1;
      for idim = 1:ndim
        Caux = kron (Proj{1,idim}(:,sub_coarse{idim}(ii)), Caux);
      end
      [ir, ic, iv] = find (Caux);
      rows(ncounter+(1:numel(ir))) = ir;
      cols(ncounter+(1:numel(ir))) = ind_coarse(ii);
      vals(ncounter+(1:numel(ir))) = iv;
      ncounter = ncounter + numel (ir);
    end
    rows = rows(1:ncounter); cols = cols(1:ncounter); vals = vals(1:ncounter);
    
    if (strcmpi (sp_coarse.space_type, 'NURBS'))
      weights_coarse = sp_coarse.weights(:);
      weights_fine = sp_fine.weights(:);
      vals = weights_fine(rows) .* vals .* weights_coarse(cols);
    end

  else
    rows = []; cols = []; vals = [];
    cumsum_ndof_coarse = sp_coarse.cumsum_ndof;  
    
    for icomp = 1:ncomp_param
      ind_comp = ind_coarse(ind_coarse>cumsum_ndof_coarse(icomp) & ...
                            ind_coarse<=cumsum_ndof_coarse(icomp+1)) - cumsum_ndof_coarse(icomp);

      spc_scalar = sp_coarse.scalar_spaces{icomp};
      spf_scalar = sp_fine.scalar_spaces{icomp};
      [rows_c, cols_c, vals_c] = matrix_change_two_levels__ (spc_scalar, spf_scalar, Proj(icomp,:), ind_comp);
      rows = [rows; rows_c]; cols = [cols; cols_c]; vals = [vals; vals_c];
    end
  end
  if (nargout == 1)
    C = sparse (rows, cols, vals, sp_fine.ndof, sp_coarse.ndof);
  end
end

if (nargout == 1)
  varargout{1} = C;
elseif (nargout == 3)
  varargout = {rows, cols, vals};
end

end
