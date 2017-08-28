% OP_DIV_V_Q_HIER: assemble the matrix B = [a(i,j)], b(i,j) = (q_i, div v_j), 
%  for hierarchical splines, exploiting the multilevel structure.
%
%   mat = op_div_v_q_hier (hspv, hspq, hmsh);
%   [rows, cols, values] = op_div_v_q_hier (hspv, hspq, hmsh);
%
% INPUT:
%
%   hspv:  object representing the vector-valued hierarchical space of trial functions (see hierarchical_space)
%   hspq:  object representing the scalar-valued hierarchical space of test functions  (see hierarchical_space)
%   hmsh:  object representing the hierarchical mesh (see hierarchical_mesh)
%
% OUTPUT:
%
%   mat:    assembled mass matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% The multilevel structure is exploited in such a way that only functions
%  of the same level of the active elements have to be computed. See also
%  op_gradu_gradv_hier for more details.
%
% Copyright (C) 2015, 2016, 2017 Eduardo M. Garau, Rafael Vazquez
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

function varargout = op_div_v_q_hier (hspv, hspq, hmsh)

  M = spalloc (hspq.ndof, hspv.ndof, 3*hspv.ndof);
  
  ndofs_v = 0;
  ndofs_q = 0;
  for ilev = 1:hmsh.nlevels
    ndofs_v = ndofs_v + hspv.ndof_per_level(ilev);
    ndofs_q = ndofs_q + hspq.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)
      x = cell(hmsh.msh_lev{ilev}.rdim,1);
      for idim = 1:hmsh.rdim
        x{idim} = reshape (hmsh.msh_lev{ilev}.geo_map(idim,:,:), hmsh.mesh_of_level(ilev).nqn, hmsh.nel_per_level(ilev));
      end
      spv_lev = sp_evaluate_element_list (hspv.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', false, 'divergence', true);
      spq_lev = sp_evaluate_element_list (hspq.space_of_level(ilev), hmsh.msh_lev{ilev});
      M_lev = op_div_v_q (spv_lev, spq_lev, hmsh.msh_lev{ilev});

      dofs_v = 1:ndofs_v;
      dofs_q = 1:ndofs_q;
      M(dofs_q,dofs_v) = M(dofs_q,dofs_v) + hspq.Csub{ilev}.' * M_lev * hspv.Csub{ilev};
    end
  end
  
  if (nargout == 1)
    varargout{1} = M;
  elseif (nargout == 3)
    [rows, cols, vals] = find (M);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end
end
