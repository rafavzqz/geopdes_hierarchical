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

function C = matrix_basis_change_new__ (hspace, lev, ind_coarse, ind_fine)

if (isa (hspace.space_of_level(1), 'sp_scalar'))
  is_scalar = true;
elseif (isa (hspace.space_of_level(1), 'sp_vector'))
  is_scalar = false;
  ncomp_param = hspace.space_of_level(lev).ncomp_param;
else
  error ('Unknown space type')
end
ndim = size (hspace.Proj{lev-1}, 2);


if (nargin < 3)
  if (is_scalar)
    C = 1;
    for idim = 1:ndim
      C = kron (hspace.Proj{lev-1,idim}, C);
    end
  else
    for icomp = 1:hspace.space_of_level(1).ncomp_param
      Caux{icomp} = 1;
      for idim = 1:ndim
        Caux{icomp} = kron (hspace.Proj{lev-1,icomp,idim}, Caux{icomp});
      end
      C = blkdiag (Caux{:});
    end
  end

elseif (nargin <= 4)
  if (is_scalar)
    sub_coarse = cell (ndim, 1);
    [sub_coarse{:}] = ind2sub ([hspace.space_of_level(lev-1).ndof_dir, 1], ind_coarse);
  
    rows = zeros (prod (hspace.space_of_level(lev).degree+1)*numel(ind_coarse), 1); cols = rows; vals = rows;
    ncounter = 0;
    for ii = 1:numel(ind_coarse)
      Caux = 1;
      for idim = 1:ndim
        Caux = kron (hspace.Proj{lev-1}{idim}(:,sub_coarse{idim}(ii)), Caux);
      end
      [ir, ic, iv] = find (Caux);
      if (nargin < 4)
          rows(ncounter+(1:numel(ir))) = ir;
          cols(ncounter+(1:numel(ir))) = ind_coarse(ii);
          vals(ncounter+(1:numel(ir))) = iv;
          ncounter = ncounter + numel (ir);
      elseif (nargin == 4)
          [~,IA,IB] = intersect(ir,ind_fine);
          rows(ncounter+(1:numel(IB))) = IB;
          cols(ncounter+(1:numel(IB))) = ii;
          vals(ncounter+(1:numel(IB))) = iv(IA);
          ncounter = ncounter + numel (IB);
      end      
    end
    rows = rows(1:ncounter); cols = cols(1:ncounter); vals = vals(1:ncounter);
  else
    rows = []; cols = []; vals = [];
    sub_coarse = cell (ndim, 1);
    cumsum_ndof_coarse = hspace.space_of_level(lev-1).cumsum_ndof;  
    cumsum_ndof_fine = hspace.space_of_level(lev).cumsum_ndof;  
    index_global = 0; counter_comp = cell(1,ncomp_param);
    for icomp = 1:ncomp_param
      ind_comp_global{icomp} = ind_coarse(ind_coarse>cumsum_ndof_coarse(icomp) & ...
                            ind_coarse<=cumsum_ndof_coarse(icomp+1));  
      ind_comp = ind_comp_global{icomp} - cumsum_ndof_coarse(icomp);
      index_global = [index_global, numel(ind_comp)];
      index_cumsum = cumsum(index_global);
      [sub_coarse{:}] = ind2sub ([hspace.space_of_level(lev-1).ndof_dir(icomp,:), 1], ind_comp);

      rows_c{icomp} = zeros (prod (hspace.space_of_level(lev).scalar_spaces{icomp}.degree+1)*numel(ind_comp), 1); vals_c{icomp} = rows_c{icomp};
      ncounter = 0;
      for ii = 1:numel(ind_comp)
        Caux = 1;
        for idim = 1:ndim
          Caux = kron (hspace.Proj{lev-1}{icomp,idim}(:,sub_coarse{idim}(ii)), Caux);
        end
        [ir, ic, iv] = find (Caux);
        if (nargin < 4)
            rows_c{icomp}(ncounter+(1:numel(ir))) = ir + cumsum_ndof_fine(icomp);
            cols_c{icomp}(ncounter+(1:numel(ir))) = ind_comp(ii) + cumsum_ndof_coarse(icomp);
            vals_c{icomp}(ncounter+(1:numel(ir))) = iv;
            ncounter = ncounter + numel (ir);
        elseif (nargin == 4)
            [~,IA,IB] = intersect(ir + cumsum_ndof_fine(icomp),ind_fine);
            counter_comp{icomp} = [counter_comp{icomp}, numel(IB)];
            rows_c{icomp}(ncounter+(1:numel(IB))) = IB;
            vals_c{icomp}(ncounter+(1:numel(IB))) = iv(IA);
            ncounter = ncounter + numel (IB);
        end
      end
%       rows = [rows; rows_c(1:ncounter)]; cols = [cols; cols_c(1:ncounter)]; vals = [vals; vals_c(1:ncounter)];  
      rows_c{icomp} = rows_c{icomp}(1:ncounter);
      vals_c{icomp} = vals_c{icomp}(1:ncounter);
    end
    i_counter = zeros(1,ncomp_param);
    ii = 1;
    for ind = ind_coarse'
        for icomp = 1:ncomp_param
            ind_c = find(ind==ind_comp_global{icomp}, 1);
            if(~isempty(ind_c))
                step = counter_comp{icomp}(ind_c);
                range = i_counter(icomp)+1:(i_counter(icomp)+step);
                cols_i(1:step,1) = ii;
                rows = [rows; rows_c{icomp}(range)]; vals = [vals; vals_c{icomp}(range)];
                i_counter(icomp) = i_counter(icomp) + step; ii = ii + 1; 
                cols = [cols; cols_i]; cols_i = [];
                break;
            end
        end
    end
  end
  if (nargin < 4)
    C = sparse (rows, cols, vals, hspace.space_of_level(lev).ndof, hspace.space_of_level(lev-1).ndof);
  elseif (nargin == 4)
    C = sparse (rows, cols, vals, numel(ind_fine), numel(ind_coarse));      
  end
end

% Computation for NURBS spaces
if (is_scalar)
  if (strcmpi (hspace.space_of_level(1).space_type, 'NURBS'))
%     Wlev = spdiags (hspace.space_of_level(lev-1).weights(:), 0, hspace.space_of_level(lev-1).ndof, hspace.space_of_level(lev-1).ndof);
%     Wlev_fine = spdiags (1./hspace.space_of_level(lev).weights(:), 0, hspace.space_of_level(lev).ndof, hspace.space_of_level(lev).ndof);
    Wlev = spdiags (hspace.space_of_level(lev-1).weights(ind_coarse), 0, numel(ind_coarse), numel(ind_coarse));
    Wlev_fine = spdiags (1./hspace.space_of_level(lev).weights(ind_fine), 0, numel(ind_fine), numel(ind_fine));
    C = Wlev_fine * C * Wlev;
  end
else
  flag = 0;
  for icomp = 1:hspace.space_of_level(1).ncomp_param
    if (strcmpi (hspace.space_of_level(lev-1).scalar_spaces{icomp}.space_type, 'NURBS'))
      flag = 1;
%       Wlev{icomp} = spdiags (hspace.space_of_level(lev-1).scalar_spaces{icomp}.weights(:), 0, ...
%           hspace.space_of_level(lev-1).scalar_spaces{icomp}.ndof, hspace.space_of_level(lev-1).scalar_spaces{icomp}.ndof);
%       Wlev_fine{icomp} = spdiags (1./hspace.space_of_level(lev).scalar_spaces{icomp}.weights(:), 0, ...
%           hspace.space_of_level(lev).scalar_spaces{icomp}.ndof, hspace.space_of_level(lev).scalar_spaces{icomp}.ndof);
      ind_comp_local_fine =  ind_fine(ind_fine>cumsum_ndof_fine(icomp) & ...
                              ind_fine<=cumsum_ndof_fine(icomp+1)) - cumsum_ndof_fine(icomp);
      ind_comp_local_coarse = ind_comp_global{icomp} - cumsum_ndof_coarse(icomp);  
      
      Wlev{icomp} = spdiags (hspace.space_of_level(lev-1).scalar_spaces{icomp}.weights(ind_comp_local_coarse), 0, ...
                    numel(ind_comp_global{icomp}), numel(ind_comp_global{icomp}));
      Wlev_fine{icomp} = spdiags (1./hspace.space_of_level(lev).scalar_spaces{icomp}.weights(ind_comp_local_fine), 0, ...
          numel(ind_comp_local_fine), numel(ind_comp_local_fine));
    else
      Wlev{icomp} = speye (hspace.space_of_level(lev-1).scalar_spaces{icomp}.ndof);
      Wlev_fine{icomp} = speye (hspace.space_of_level(lev).scalar_spaces{icomp}.ndof);
    end
  end
  if (flag)
    Wlev = blkdiag (Wlev{:});
    Wlev_fine = blkdiag (Wlev_fine{:});
    C = Wlev_fine * C * Wlev;
  end
end

if (hspace.truncated)
  indices = union (hspace.active{lev}, hspace.deactivated{lev});
  C(indices,:) = 0;
end

end
