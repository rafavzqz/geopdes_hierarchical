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

npatch = hspace.space_of_level(1).npatch;

sp_coarse = hspace.space_of_level(lev-1);
sp_fine = hspace.space_of_level(lev);
Proj = hspace.Proj(lev-1,:);

for iptc = 1:npatch
  spc_patch = sp_coarse.sp_patch{iptc};
  spf_patch = sp_fine.sp_patch{iptc};
  matrix_change_two_levels__ (spc_patch, spf_patch
end


if (isa (hspace.space_of_level(1).sp_patch{1}, 'sp_scalar'))
  is_scalar = true;
  ndim = size (hspace.Proj{1}, 2);
elseif (isa (hspace.space_of_level(1).sp_patch{1}, 'sp_vector'))
  is_scalar = false;
  ndim = size (hspace.Proj{1}, 2);
else
  error ('Unknown space type')
end

npatch = hspace.space_of_level(1).npatch;
ndim = size (hspace.Proj{1}, 2);

if (nargin < 3)
  C = sparse (hspace.space_of_level(lev).ndof, hspace.space_of_level(lev-1).ndof);
  if (is_scalar) % Full matrix, scalar case
    for ipatch = 1:npatch
      Cpatch = 1;
      Proj = hspace.Proj{lev-1, ipatch};
      for idim = 1:ndim
        Cpatch = kron (Proj{idim}, Cpatch);
      end

      if (strcmpi (hspace.space_of_level(1).sp_patch{ipatch}.space_type, 'NURBS'))
        sp_patch = hspace.space_of_level(lev-1).sp_patch{ipatch};
        Wlev = spdiags (sp_patch.weights(:), 0, sp_patch.ndof, sp_patch.ndof);
        sp_patch = hspace.space_of_level(lev).sp_patch{ipatch};
        Wlev_fine = spdiags (1./sp_patch.weights(:), 0, sp_patch.ndof, sp_patch.ndof);
        Cpatch = Wlev_fine * Cpatch * Wlev;
      end
      C(hspace.space_of_level(lev).gnum{ipatch}, hspace.space_of_level(lev-1).gnum{ipatch}) = Cpatch;
    end
  else % Full matrix, vector case
    for ipatch = 1:npatch
      for icomp = 1:hspace.space_of_level(1).sp_patch{ipatch}.ncomp_param
        Caux{icomp} = 1;
        Proj = hspace.Proj{lev-1, ipatch};
        for idim = 1:ndim
          Caux{icomp} = kron (Proj{icomp,idim}, Caux{icomp});
        end
        Cpatch = blkdiag (Caux{:});
      end
% XXXXXXXXXXXXXXXXXXXX APPLY THE NURBS STUFF
      C(hspace.space_of_level(lev).gnum{ipatch}, hspace.space_of_level(lev-1).gnum{ipatch}) = Cpatch;
    end
      
  end

elseif (nargin == 3)
  nrows = hspace.space_of_level(lev).ndof; ncols = hspace.space_of_level(lev-1).ndof;
  if (is_scalar) % Not full matrix, scalar version
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

      if (strcmpi (hspace.space_of_level(1).sp_patch{ipatch}.space_type, 'NURBS'))
        sp_patch = hspace.space_of_level(lev-1).sp_patch{ipatch};
        Wcoarse = sp_patch.weights(:);
        sp_patch = hspace.space_of_level(lev).sp_patch{ipatch};
        Wfine = 1./sp_patch.weights(:);
        vals(nc_in_patch:ncounter) = vals(nc_in_patch:ncounter) .* Wcoarse(rows(nc_in_patch:ncounter)) .* Wfine(rows(nc_in_patch:ncounter));
      end
      rows(nc_in_patch:ncounter) = hspace.space_of_level(lev).gnum{ipatch}(rows(nc_in_patch:ncounter));
      cols(nc_in_patch:ncounter) = hspace.space_of_level(lev-1).gnum{ipatch}(cols(nc_in_patch:ncounter));
    end
  else % Not full matrix, vector version
    ncomp_param = hspace.space_of_level(1).sp_patch{1}.ncomp_param;
    size_alloc = numel (ind_coarse) * prod (hspace.space_of_level(lev).sp_patch{1}.scalar_spaces{1}.degree + 1) * ncomp_param;
    rows = zeros (size_alloc, 1); cols = rows; vals = rows;
    ncounter = 0;
    for ipatch = 1:npatch
      gnum_coarse = hspace.space_of_level(lev-1).gnum{ipatch};
      [ind_ptc,~,local_indices] = intersect (ind_coarse, gnum_coarse);
      gnum_fine = hspace.space_of_level(lev).gnum{ipatch};
    
      sub_coarse = cell (ndim, 1);
      cumsum_ndof_coarse = hspace.space_of_level(lev-1).sp_patch{ipatch}.cumsum_ndof;  
      cumsum_ndof_fine = hspace.space_of_level(lev).sp_patch{ipatch}.cumsum_ndof;
      for icomp = 1:ncomp_param
        nc_in_comp = ncounter+1;
        scalar_space = hspace.space_of_level(lev).sp_patch{ipatch}.scalar_spaces{icomp};
        ind_comp = local_indices(local_indices>cumsum_ndof_coarse(icomp) & ...
                                 local_indices<=cumsum_ndof_coarse(icomp+1)) - cumsum_ndof_coarse(icomp);
        [sub_coarse{:}] = ind2sub ([hspace.space_of_level(lev-1).sp_patch{ipatch}.ndof_dir(icomp,:), 1], ind_comp);

        rows_c = zeros (prod (scalar_space.degree+1)*numel(ind_comp), 1); cols_c = rows_c; vals_c = rows_c;
        nc_comp = 0;
        for ii = 1:numel(ind_comp)
          Proj = hspace.Proj{lev-1, ipatch};
          Caux = 1;
          for idim = 1:ndim
            Caux = kron (Proj{icomp,idim}(:,sub_coarse{idim}(ii)), Caux);
 %             Caux = kron (hspace.Proj{lev-1,icomp,idim}(:,sub_coarse{idim}(ii)), Caux);
          end
          [ir, ic, iv] = find (Caux);
          rows_c(nc_comp+(1:numel(ir))) = ir + cumsum_ndof_fine(icomp);
          cols_c(nc_comp+(1:numel(ir))) = ind_comp(ii) + cumsum_ndof_coarse(icomp);
          vals_c(nc_comp+(1:numel(ir))) = iv;
          nc_comp = nc_comp + numel (ir);
          ncounter = ncounter + numel (ir);
        end
        
        rows(nc_in_comp:ncounter) = gnum_fine(rows_c(1:nc_comp));
        cols(nc_in_comp:ncounter) = gnum_coarse(cols_c(1:nc_comp));
        vals(nc_in_comp:ncounter) = vals_c(1:nc_comp);
% XXXXXXXXXXXXXXXXXXXX APPLY THE NURBS STUFF
      end
    end

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
