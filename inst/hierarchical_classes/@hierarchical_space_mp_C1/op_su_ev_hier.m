% OP_SU_EV_HIER: assemble the matrix A = [a(i,j)], a(i,j) = 1/2 (sigma (u_j), epsilon (v_i))
%  for hierarchical splines, exploiting the multilevel structure.
%
%   mat = op_su_ev_hier (hspu, hspv, hmsh, lambda, mu);
%   [rows, cols, values] = op_su_ev_hier (hspu, hspv, hmsh, lambda, mu);
%
% INPUT:
%
%   hspu:  object representing the hierarchical space of trial functions (see hierarchical_space_mp_C1)
%   hspv:  object representing the hierarchical space of test functions  (see hierarchical_space_mp_C1)
%   hmsh:  object representing the hierarchical mesh (see hierarchical_mesh_mp)
%   lambda, mu: function handles to compute the Lame' coefficients
%   patch_list: list of patches on which to compute the integrals
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
% Copyright (C) 2022, 2023 Rafael Vazquez
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

function varargout = op_su_ev_hier (hspu, hspv, hmsh, lambda, mu, patch_list)

  if (nargin < 6)
    patch_list = 1:hmsh.mesh_of_level(1).npatch;
  end
  
  rdim = hmsh.rdim;
  K = spalloc (rdim*hspv.ndof, rdim*hspu.ndof, 3*hspu.ndof);

  ndofs_u = 0;
  ndofs_v = 0;
  for ilev = 1:hmsh.nlevels
    ndofs_u = ndofs_u + hspu.ndof_per_level(ilev);
    ndofs_v = ndofs_v + hspv.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)
      msh_lev = msh_restrict_to_patches (hmsh.msh_lev{ilev}, patch_list);
        
      if (msh_lev.nel > 0)
        x = cell (rdim,1);
        for idim = 1:rdim
          x{idim} = reshape (msh_lev.geo_map(idim,:,:), msh_lev.nqn, msh_lev.nel);
        end
        spu_lev = sp_evaluate_element_list (hspu.space_of_level(ilev), msh_lev, 'value', false, 'gradient', true);
        spv_lev = sp_evaluate_element_list (hspv.space_of_level(ilev), msh_lev, 'value', false, 'gradient', true);

        spu_lev = change_connectivity_localized_Csub (spu_lev, hspu, ilev);
        spv_lev = change_connectivity_localized_Csub (spv_lev, hspv, ilev);
        
        spu_lev = sp_scalar2vector (spu_lev, msh_lev, 'value', false, 'gradient', true, 'divergence', true);
        spv_lev = sp_scalar2vector (spv_lev, msh_lev, 'value', false, 'gradient', true, 'divergence', true);
        K_lev = op_su_ev (spu_lev, spv_lev, msh_lev, lambda (x{:}), mu (x{:}));

        dofs_u = [];
        dofs_v = [];
        for icomp = 1:rdim
          dofs_u = union (dofs_u, (icomp-1)*hspu.ndof + (1:ndofs_u));
          dofs_v = union (dofs_v, (icomp-1)*hspv.ndof + (1:ndofs_v));
        end

        Csub_u = repmat (hspu.Csub(ilev), 1, rdim);
        Csub_u = blkdiag (Csub_u{:});
        Csub_v = repmat (hspv.Csub(ilev), 1, rdim);
        Csub_v = blkdiag (Csub_v{:});

        K(dofs_v,dofs_u) = K(dofs_v,dofs_u) + Csub_v.' * K_lev * Csub_u;
      end
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
