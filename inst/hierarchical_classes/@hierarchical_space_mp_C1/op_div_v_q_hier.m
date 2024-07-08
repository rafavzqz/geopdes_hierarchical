% OP_DIV_V_Q_HIER: assemble the matrix B = [a(i,j)], b(i,j) = (q_i, div v_j), 
%  for hierarchical splines, exploiting the multilevel structure.
%
%   mat = op_div_v_q_hier (hspv, hspq, hmsh);
%   [rows, cols, values] = op_div_v_q_hier (hspv, hspq, hmsh);
%
% INPUT:
%
%   hspv:  object representing the vector-valued hierarchical space of trial functions (see hierarchical_space_mp_C1)
%   hspq:  object representing the scalar-valued hierarchical space of test functions  (see hierarchical_space_mp_C1)
%   hmsh:  object representing the hierarchical mesh (see hierarchical_mesh_mp)
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
% Copyright (C) 2015, 2016, 2018 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2022 Rafael Vazquez
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

function varargout = op_div_v_q_hier (hspv, hspq, hmsh, patch_list)

  if (nargin < 4)
    patch_list = 1:hmsh.mesh_of_level(1).npatch;
  end
  
  M = spalloc (hspq.ndof, hspv.ndof, 3*hspv.ndof);
  
  ndofs_v = 0;
  ndofs_q = 0;
  for ilev = 1:hmsh.nlevels
    ndofs_v = ndofs_v + hspv.ndof_per_level(ilev);
    ndofs_q = ndofs_q + hspq.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)
      msh_lev = msh_restrict_to_patches (hmsh.msh_lev{ilev}, patch_list);

      if (msh_lev.nel > 0)
        x = cell(msh_lev.rdim,1);
        for idim = 1:hmsh.rdim
          x{idim} = reshape (msh_lev.geo_map(idim,:,:), msh_lev.nqn, msh_lev.nel);
        end
        spv_lev = sp_evaluate_element_list (hspv.space_of_level(ilev), msh_lev, 'value', false, 'divergence', true);
        spq_lev = sp_evaluate_element_list (hspq.space_of_level(ilev), msh_lev);

        spv_lev = change_connectivity_localized_Csub (spv_lev, hspv, ilev);
        spq_lev = change_connectivity_localized_Csub (spq_lev, hspq, ilev);
        M_lev = op_u_v (spv_lev, spq_lev, msh_lev, coeff (x{:}));

        dofs_v = 1:ndofs_v;
        dofs_q = 1:ndofs_q;
        M(dofs_q,dofs_v) = M(dofs_q,dofs_v) + hspq.Csub{ilev}.' * M_lev * hspv.Csub{ilev};
      end
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
