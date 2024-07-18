% OP_PENALTY_DUDN: matrix and right-hand side to compute penalty terms to impose du/dn = f.
%  It computes the terms of the form Cp*(du/dn, dv/dn) and Cp*(f, dv/dn), 
%  using the same space for trial and test functions.
%
%   [mat, rhs] = op_penalty_dudn (hspace, hmsh, sides, Cpen, [coeff]);
%
% INPUT:
%
%  hspace: object representing the space of trial functions (see hierarchical_space)
%  hmsh:   object defining the domain partition and the quadrature rule (see hierarchical_mesh)
%  sides:  boundary sides on which to compute the integrals
%  Cpen:   penalization parameter, Cp in the equation above
%  coeff:  function handle to compute the Neumann condition. If empty, the returned rhs will be zero.
%
% OUTPUT:
%
%  mat:    assembled mass matrix
%  rhs:    assembled right-hand side
% 
% Copyright (C) 2023, 2024 Michele Torre, Rafael Vazquez
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

function [mass_pen, rhs] = op_penalty_dudn (hspace, hmsh, nmnn_sides, Cpen, coeff)
    
  mass_pen =  spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);
  rhs = zeros(hspace.ndof,1);
  
  for iside = nmnn_sides
    hmsh_side_int = hmsh_boundary_side_from_interior (hmsh, iside);
    shifting_indices = cumsum ([0 hmsh.boundary(iside).nel_per_level]);
  
    ndofs_u = 0;
  
    for ilev = 1:hmsh.boundary(iside).nlevels
      ndofs_u = ndofs_u + hspace.ndof_per_level(ilev);
  
      if (hmsh.boundary(iside).nel_per_level(ilev) > 0)
  
        % mesh of the selected side
        elements = shifting_indices(ilev)+1:shifting_indices(ilev+1);
        hmsh_side = hmsh_eval_boundary_side (hmsh, iside, elements);
  
        % reconstruct the part of the space defined on the selected boundary
        msh_side_from_interior = hmsh_side_int.mesh_of_level(ilev);
        sp_bnd = hspace.space_of_level(ilev).constructor (msh_side_from_interior);
        msh_side_from_interior_struct = msh_evaluate_element_list (msh_side_from_interior, hmsh_side_int.active{ilev});
  
        % evaluate the space
        spu_lev = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'gradient', true);
        spu_lev = change_connectivity_localized_Csub (spu_lev, hspace, ilev);
  
        % compute the matrix
        coe_side = repmat(hmsh_side.element_size, hmsh_side.nqn,1);
        coe_side = Cpen .* coe_side;
        K_lev = op_gradu_n_gradv_n(spu_lev, spu_lev, hmsh_side, coe_side);
  
        dofs_u = 1:ndofs_u;
        Ktmp =  hspace.Csub{ilev}.' * K_lev * hspace.Csub{ilev};
        mass_pen(dofs_u,dofs_u) = mass_pen(dofs_u,dofs_u) + Ktmp;

        if (nargin == 5)
          x = cell (hmsh.rdim, 1);
          for idim = 1:hmsh.rdim
            x{idim} = reshape (hmsh_side.geo_map(idim,:,:), hmsh_side.nqn, hmsh_side.nel);
          end
          rhs_lev = op_gradv_n_f (spu_lev, hmsh_side, Cpen*coeff(x{:}));
          rhs(dofs_u) = rhs(dofs_u) + hspace.Csub{ilev}.' * rhs_lev;
        end
      end
    end
  end

end
