% OP_KL_SHELLS_HIER: assemble the Kirchhoff-Love shell stiffness matrix.
%
%   mat = op_KL_shells_hier (hspu, hspv, hmsh, E_coeff, nu_coeff, t_coeff);
%
% INPUT:
%
%  hspu:     object representing the space of trial functions (see hierarchical_space_mp_C1)
%  hspu:     object representing the space of test functions (see hierarchical_space_mp_C1)
%  hmsh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%  E_coeff:  function handle to compute the Young's modulus
%  nu_coeff: function handle to compute the Poisson's ratio
%  t_coeff:  thickness of the shell, scalar value
%
% OUTPUT:
%
%  mat:    assembled stiffness matrix
% 
% Copyright (C) 2023 Rafael Vazquez
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

function varargout = op_KL_shells_hier (hspu, hspv, hmsh, E_coeff, nu_coeff, t_coeff)

  K = sparse (hspv.ndof, hspu.ndof);

  ndofs_u = 0;
  ndofs_v = 0;
  for ilev = 1:hmsh.nlevels
    ndofs_u = ndofs_u + hspu.ndof_per_level(ilev);
    ndofs_v = ndofs_v + hspv.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)
      x = cell (hmsh.rdim,1);
      for idim = 1:hmsh.rdim
        x{idim} = reshape (hmsh.msh_lev{ilev}.geo_map(idim,:,:), hmsh.mesh_of_level(ilev).nqn, hmsh.nel_per_level(ilev));
      end
      spu_lev = sp_evaluate_element_list_param (hspu.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true, 'gradient', true, 'hessian', true);
      spv_lev = sp_evaluate_element_list_param (hspv.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true, 'gradient', true, 'hessian', true);

      spu_lev = change_connectivity_localized_Csub (spu_lev, hspu, ilev);
      spv_lev = change_connectivity_localized_Csub (spv_lev, hspv, ilev);

      K_lev = op_KL_shells (spu_lev, spv_lev, hmsh.msh_lev{ilev}, E_coeff(x{:}), nu_coeff(x{:}), t_coeff);

      dofs_u = 1:ndofs_u;
      dofs_v = 1:ndofs_v;
      K(dofs_v,dofs_u) = K(dofs_v,dofs_u) + hspv.Csub{ilev}.' * K_lev * hspu.Csub{ilev};
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
