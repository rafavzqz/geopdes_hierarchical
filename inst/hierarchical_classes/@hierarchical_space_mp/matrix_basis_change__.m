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
% Copyright (C) 2017-2019 Rafael Vazquez
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

C = sparse (sp_fine.ndof, sp_coarse.ndof);

if (nargin == 3)
  for iptc = 1:npatch
    spc_patch = sp_coarse.sp_patch{iptc};
    spf_patch = sp_fine.sp_patch{iptc};
    [~,local_indices,~] = intersect (sp_coarse.gnum{iptc}, ind_coarse);
    Cpatch = matrix_change_two_levels__ (spc_patch, spf_patch, Proj{iptc}, local_indices);
    C(sp_fine.gnum{iptc},sp_coarse.gnum{iptc}) = Cpatch;
  end
else
  for iptc = 1:npatch
    spc_patch = sp_coarse.sp_patch{iptc};
    spf_patch = sp_fine.sp_patch{iptc};
    Cpatch = matrix_change_two_levels__ (spc_patch, spf_patch, Proj{iptc});
    C(sp_fine.gnum{iptc},sp_coarse.gnum{iptc}) = Cpatch;
  end
end

if (hspace.truncated)
  indices = union (hspace.active{lev}, hspace.deactivated{lev});
  C(indices,:) = 0;
end

end
