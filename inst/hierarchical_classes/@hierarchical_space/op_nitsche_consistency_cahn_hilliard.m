% OP_NITSCHE_CONSISTENCY_CAHN_HILLIARD: assemble the matrix A = [a(i,j)], a(i,j) = (epsilon (grad v n)_j, Delta u_i), 
%  with n the normal vector to the boundary, using the same space for trial and test functions.
%
%   mat = op_nitsche_consistency_cahn_hilliard (hspace, hmsh, sides, [coeff]);
%
% INPUT:
%
%  hspace: object representing the space of trial functions (see hierarchical_space)
%  hmsh:   object defining the hierarchical mesh (see hierarchical_mesh)
%  sides:  boundary sides on which to compute the integrals
%  coeff:  function handle to compute the epsilon coefficient. If empty, it is taken equal to one.
%
% OUTPUT:
%
%  mat:    assembled matrix
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

function A = op_nitsche_consistency_cahn_hilliard (hspace, hmsh, sides, coeff)

  if (nargin == 3 || isempty(coeff))
    coeff = @(varargin) ones (size(varargin{1}));
  end

  if (~isempty (sides))
    A =  spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);
    for iside = sides
      hmsh_side_int = hmsh_boundary_side_from_interior (hmsh, iside);
      shifting_indices = cumsum ([0 hmsh.boundary(iside).nel_per_level]);

      ndofs = 0;
      for ilev = 1:hmsh.boundary(iside).nlevels

        ndofs = ndofs + hspace.ndof_per_level(ilev);

        if (hmsh.boundary(iside).nel_per_level(ilev) > 0)
          % mesh of the selected side
          elements = shifting_indices(ilev)+1:shifting_indices(ilev+1);
          hmsh_side = hmsh_eval_boundary_side (hmsh, iside, elements);

          x = cell (hmsh.rdim,1);
          for idim = 1:hmsh.rdim
            x{idim} = reshape (hmsh_side.geo_map(idim,:,:), hmsh_side.nqn, hmsh_side.nel);
          end

          % reconstruct the part of the space defined on the selected boundary
          msh_side_from_interior = hmsh_side_int.mesh_of_level(ilev);
          sp_bnd = hspace.space_of_level(ilev).constructor (msh_side_from_interior);
          msh_side_from_interior_struct = msh_evaluate_element_list (msh_side_from_interior, hmsh_side_int.active{ilev});

          % evaluate the space
          sp_lev = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', false, 'gradient', true, 'laplacian', true);
          sp_lev = change_connectivity_localized_Csub (sp_lev, hspace, ilev);

          % compute the matrix
          K_lev = op_gradv_n_laplaceu (sp_lev, sp_lev, hmsh_side, coeff(x{:}));

          dofs = 1:ndofs;
          A(dofs,dofs) = A(dofs,dofs) + hspace.Csub{ilev}.' * K_lev * hspace.Csub{ilev};
        end
      end

    end
  else
    A = [];
  end

end
