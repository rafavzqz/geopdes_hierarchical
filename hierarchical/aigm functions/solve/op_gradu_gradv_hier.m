% OP_GRADU_GRADV_HIER: assemble the stiffness matrix A = [a(i,j)], a(i,j) = (epsilon grad u_j, grad v_i), 
%  for hierarchical splines, exploiting the multilevel structure.
%
%   mat = op_gradu_gradv_hier (hspu, hspv, hmsh, epsilon);
%   [rows, cols, values] = op_gradu_gradv_hier (hspu, hspv, hmsh, epsilon);
%
% INPUT:
%
%   hspu:    object representing the hierarchical space of trial functions (see XXXXXXXXX)
%   hspv:    object representing the hierarchical space of test functions  (see XXXXXXXXX)
%   hmsh:    object representing the hierarchical mesh (see XXXXXXXXXXXXX)
%   epsilon: function handle to compute the diffusion coefficient
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% For each level, a matrix is computed from the integrals on the active
%  elements of that level, using only functions of that level (both active
%  and inactive). Then, the matrix for the hierarchical space is computed
%  using the matrix relation between the basis functions of different levels.
%
% Copyright (C) 2015, Eduardo M. Garau, Rafael Vazquez
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

function varargout = op_gradu_gradv_hier (hspu, hspv, hmsh, coeff)

  K = spalloc (hspv.ndof, hspu.ndof, 3*hspu.ndof);

  ndof_per_level = hspu.ndof_per_level;
  dif = hmsh.nlevels - hspu.nlevels;
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%  CAN WE REMOVE THE WARNING (and dif)?
  if (dif)
    warning('Different number of levels for mesh and space. This should never happen')
    ndof_per_level = [ndof_per_level(:); zeros(dif,1)];
  end

  ndofs_u = 0;
  ndofs_v = 0;
  for ilev = 1:hmsh.nlevels
    ndofs_u = ndofs_u + ndof_per_level(ilev);%hspu.ndof_per_level(ilev);
    ndofs_v = ndofs_v + ndof_per_level(ilev);%hspv.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)
      x = cell (hmsh.rdim,1);
      for idim = 1:hmsh.rdim
        x{idim} = hmsh.msh_lev{ilev}.geo_map(idim,:,:);
      end
      K_lev = op_gradu_gradv (hspu.sp_lev{ilev}, hspv.sp_lev{ilev}, hmsh.msh_lev{ilev}, coeff (x{:}));

      dofs_u = 1:ndofs_u;
      dofs_v = 1:ndofs_v;

      K(dofs_v,dofs_u) = K(dofs_v,dofs_u) + hspv.C{ilev}'*K_lev*hspu.C{ilev};
    end
  end

  if (nargout == 1)
    varargout{1} = K;
  elseif (nargout == 3)
    [rows, cols, vals] = find (K);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end
