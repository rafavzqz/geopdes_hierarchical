% OP_EU_EV_HIER: assemble the matrix A = [a(i,j)], a(i,j) = 2*mu (epsilon (u_j), epsilon (v_i))
%  for hierarchical splines, exploiting the multilevel structure.
%
%   mat = op_eu_ev_hier (hspu, hspv, hmsh, mu);
%   [rows, cols, values] = op_eu_ev_hier (hspu, hspv, hmsh, mu);
%
% INPUT:
%
%   hspu:  object representing the hierarchical space of trial functions (see hierarchical_space)
%   hspv:  object representing the hierarchical space of test functions  (see hierarchical_space)
%   hmsh:  object representing the hierarchical mesh (see hierarchical_mesh)
%   mu:    function handle to compute the Lame' coefficient
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% The multilevel structure is exploited in such a way that only functions
%  of the same level of the active elements have to be computed. See also
%  op_gradu_gradv_hier for more details.
%
% Copyright (C) 2015, 2016, 2018 Eduardo M. Garau, Rafael Vazquez
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

function varargout = op_eu_ev_hier (hspu, hspv, hmsh, mu, patch_list)

  if (nargin < 5)
    patch_list = 1:hmsh.mesh_of_level(1).npatch;
  end
  
  M = spalloc (hspv.ndof, hspu.ndof, 3*hspu.ndof);
  
  ndofs_u = 0;
  ndofs_v = 0;
  for ilev = 1:hmsh.nlevels
    ndofs_u = ndofs_u + hspu.ndof_per_level(ilev);
    ndofs_v = ndofs_v + hspv.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)
      msh_lev = msh_restrict_to_patches (hmsh.msh_lev{ilev}, patch_list);
        
      if (msh_lev.nel > 0)
        x = cell(msh_lev.rdim,1);
        for idim = 1:hmsh.rdim
          x{idim} = reshape (msh_lev.geo_map(idim,:,:), msh_lev.nqn, msh_lev.nel);
        end
        spu_lev = sp_evaluate_element_list (hspu.space_of_level(ilev), msh_lev, 'value', false, 'gradient', true);
        spv_lev = sp_evaluate_element_list (hspv.space_of_level(ilev), msh_lev, 'value', false, 'gradient', true);
        M_lev = op_eu_ev (spu_lev, spv_lev, hmsh.msh_lev{ilev}, mu (x{:}));

        dofs_u = 1:ndofs_u;
        dofs_v = 1:ndofs_v;
        M(dofs_v,dofs_u) = M(dofs_v,dofs_u) + hspv.Csub{ilev}.' * M_lev * hspu.Csub{ilev};
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
