% % % % % % OP_GRADU_GRADV_HIER: assemble the stiffness matrix A = [a(i,j)], a(i,j) = (epsilon grad u_j, grad v_i), 
% % % % % %  for hierarchical splines, exploiting the multilevel structure.
% % % % % %
% % % % % %   mat = op_gradu_gradv_hier (hspu, hspv, hmsh, [epsilon]);
% % % % % %   [rows, cols, values] = op_gradu_gradv_hier (hspu, hspv, hmsh, [epsilon]);
% % % % % %
% % % % % % INPUT:
% % % % % %
% % % % % %   hspu:    object representing the hierarchical space of trial functions (see hierarchical_space)
% % % % % %   hspv:    object representing the hierarchical space of test functions  (see hierarchical_space)
% % % % % %   hmsh:    object representing the hierarchical mesh (see hierarchical_mesh)
% % % % % %   epsilon: function handle to compute the diffusion coefficient
% % % % % %
% % % % % % OUTPUT:
% % % % % %
% % % % % %   mat:    assembled stiffness matrix
% % % % % %   rows:   row indices of the nonzero entries
% % % % % %   cols:   column indices of the nonzero entries
% % % % % %   values: values of the nonzero entries
% % % % % % 
% % % % % % For each level, a matrix is computed from the integrals on the active
% % % % % %  elements of that level, using only functions of that level (both active
% % % % % %  and inactive). Then, the matrix for the hierarchical space is computed
% % % % % %  using the matrix relation between the basis functions of different levels.
% % % % % %
% % % % % % Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
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

function varargout = op_v_gradp_hier (hspv, hspp, hmsh, coeff)

  if (nargin == 3)
    coeff = @(varargin) ones (size(varargin{1}));
  end

  K = spalloc (hspp.ndof, hspv.ndof, 3*hspv.ndof);

  ndofs_v = 0;
  ndofs_p = 0;
  for ilev = 1:hmsh.nlevels
    ndofs_p = ndofs_p + hspp.ndof_per_level(ilev);
    ndofs_v = ndofs_v + hspv.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)
      x = cell (hmsh.rdim,1);
      for idim = 1:hmsh.rdim
        x{idim} = reshape (hmsh.msh_lev{ilev}.geo_map(idim,:,:), hmsh.mesh_of_level(ilev).nqn, hmsh.nel_per_level(ilev));
      end
      spv_lev = sp_evaluate_element_list (hspv.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
      spp_lev = sp_evaluate_element_list (hspp.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', false, 'gradient', true);
      K_lev = op_v_gradp (spv_lev, spp_lev, hmsh.msh_lev{ilev}, coeff (x{:}));

      dofs_p = 1:ndofs_p;
      dofs_v = 1:ndofs_v;

      K(dofs_p,dofs_v) = K(dofs_p,dofs_v) + hspp.Csub{ilev}.' * K_lev * hspv.Csub{ilev};
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
