% HSPACE_SUBDIVISION_MATRIX: compute the matrices for changing basis, from
%                 active functions to B-splines (C^1 basis) of the tensor product spaces.
%
%   Csub = hspace_subdivision_matrix (hspace, [hmsh])
%
% If the hierarchical mesh is present, only the rows relative to functions
%  on active and deactivated elements are computed. The matrix size is not affected.
%
% INPUT:
%
%   hspace:    object representing the hierarchical space (see hierarchical_space_mp)
%   hmsh:      object representing the hierarchical mesh, only needed for the 'reduced' version (see hierarchical_mesh_mp)
%
% OUTPUT:
%
%   Csub:      cell-array with the matrices for basis change. The size of the matrix Csub{lev} is
%                 hspace.space_of_level(lev).ndof  x  sum(hspace.ndof_per_level(1:lev))
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2018, 2019 Luca Coradello, Rafael Vazquez
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

function [Csub, row_indices] = hspace_subdivision_matrix (hspace, hmsh, option)

if (nargin == 1)
  option = 'full';
elseif (nargin == 2)
  option = 'reduced';
elseif (~strcmpi (option, 'reduced') && ~strcmpi (option, 'full'))
  warning ('Unknown option. Trying to compute the reduced version')
end

Csub = cell (hspace.nlevels, 1);
row_indices = cell (hspace.nlevels, 1);
Csub{1} = speye (hspace.space_of_level(1).ndof); 

if (strcmpi (option, 'reduced'))
  Csub{1} = Csub{1}(:,hspace.active{1});
  fun_on_active = sp_get_basis_functions (hspace.space_of_level(1), hmsh.mesh_of_level(1), hmsh.active{1});
  fun_on_deact = sp_get_basis_functions (hspace.space_of_level(1), hmsh.mesh_of_level(1), hmsh.deactivated{1});
  fun_on_deact = union (fun_on_active, fun_on_deact);
  row_indices{1} = fun_on_deact;

  for lev = 2:hmsh.nlevels
    fun_on_active_finer_level = sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), hmsh.active{lev});
    fun_on_deact_finer_level = sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), hmsh.deactivated{lev});
    fun_on_deact_finer_level = union (fun_on_active_finer_level, fun_on_deact_finer_level);

    [~,~,IB] = intersect(hspace.active{lev}, fun_on_deact_finer_level);
    I_rows = IB; I_cols = 1:hspace.ndof_per_level(lev);
    I = sparse (I_rows, I_cols, ones(size(I_cols)), numel(fun_on_deact_finer_level), hspace.ndof_per_level(lev));

    aux = matrix_basis_change__ (hspace, lev, fun_on_deact, fun_on_deact_finer_level);

    fun_on_deact = fun_on_deact_finer_level;
    Csub{lev} = [aux*Csub{lev-1}, I];
    row_indices{lev} = fun_on_deact;
  end
  
elseif (strcmpi (option, 'full'))
  Csub{1} = Csub{1}(:,hspace.active{1});
  for lev = 2:hspace.nlevels
    I = speye (hspace.space_of_level(lev).ndof); 
    aux = matrix_basis_change__ (hspace, lev);
    Csub{lev} = [aux*Csub{lev-1}, I(:,hspace.active{lev})];
    clear aux I
  end
end

end
