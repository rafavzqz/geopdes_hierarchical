% SP_NITSCHE_KL_ROTATION: impose zero rotation for Kirchhoff-Love shells using Nitsche method
%
%   N_mat = sp_nitsche_KL_rotation (space, msh, bnd_sides, E_coeff, nu_coeff, thickness, C_penalty)
%
% INPUT:
%     
%    space:     hierarchical space object (see hierarchical_space_mp_C1)
%    msh:       hierarchical mesh object (see hierarchical_mesh_mp)
%    bnd_sides: boundary sides on which the rotation free condition is imposed
%    E_coeff:   function handle for the Young modulus
%    nu_coeff:  function handle for the Poisson ratio
%    thickness: scalar value for the thickness
%    C_penalty: parameter for the penalization term
%   
% OUTPUT:
%
%     N_mat:       the computed matrix, to be added in the left hand-side
%
% Copyright (C) 2023, 2024 Giuliano Guarino, Rafael Vazquez
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
function A = sp_nitsche_KL_rotation (hspace, hmsh, bnd_sides, E_coeff, nu_coeff, thickness, penalty_coeff)

  rdim = hmsh.rdim;
  A = spalloc (rdim*hspace.ndof, rdim*hspace.ndof, 3*rdim*hspace.ndof);

  boundaries = hmsh.mesh_of_level(1).boundaries;
  Nbnd = cumsum ([0, boundaries.nsides]);
  last_dof = cumsum (hspace.ndof_per_level);
  
% Compute the matrices to impose the tangential boundary condition weakly
  for iref = bnd_sides
    iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    ndofs_u = 0;
    for ilev = 1:hmsh.boundary.nlevels
      ndofs_u = ndofs_u + hspace.ndof_per_level(ilev);
      patch_bnd_shifting = cumsum ([0 hmsh.boundary.mesh_of_level(ilev).nel_per_patch]);
      % patch_shifting = cumsum ([0 hmsh.mesh_of_level(ilev).nel_per_patch]);
      % dofs_on_lev = 1:last_dof(ilev);

      for ii = 1:numel(iref_patch_list)
        iptc_bnd = iref_patch_list(ii);
        iptc = boundaries(iref).patches(ii);
        iside = boundaries(iref).faces(ii);
        elems_bnd_patch = patch_bnd_shifting(iptc_bnd)+1:patch_bnd_shifting(iptc_bnd+1);
        [~, ~, elements] = intersect (hmsh.boundary.active{ilev}, elems_bnd_patch);

        if (~isempty (elements))
          msh_patch = hmsh.mesh_of_level(ilev).msh_patch{iptc};

          msh_side = msh_eval_boundary_side (msh_patch, iside, elements);
          msh_side_from_interior = msh_boundary_side_from_interior (msh_patch, iside);

          sp_bnd = hspace.space_of_level(ilev).sp_patch{iptc}.constructor (msh_side_from_interior);
          msh_side_fi = msh_evaluate_element_list (msh_side_from_interior, elements);
          sp_bnd_param = sp_evaluate_element_list_param (sp_bnd, msh_side_fi, 'value', true, 'gradient', true, 'hessian', true);
          sp_bnd_param = sp_scalar2vector_param (sp_bnd_param, msh_side_fi, 'value', true, 'gradient', true, 'hessian', true);

          for idim = 1:rdim
            x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
          end
          
% Evaluating parameters on the boundary
          E_bnd  = reshape(E_coeff  (x{:}), msh_side.nqn, msh_side.nel);
          nu_bnd = reshape(nu_coeff (x{:}), msh_side.nqn, msh_side.nel);

          pen_coeff = penalty_coeff * 2^(-ilev);
          A_side = op_nitsche_KL_boundary (sp_bnd_param, sp_bnd_param, msh_side, msh_side_fi, E_bnd, nu_bnd, thickness, pen_coeff);

          [Cpatch, Cpatch_cols_lev] = sp_compute_Cpatch (hspace.space_of_level(ilev), iptc);
          [~,Csub_rows,Cpatch_cols] = intersect (hspace.Csub_row_indices{ilev}, Cpatch_cols_lev);
          Caux = Cpatch(:,Cpatch_cols) * hspace.Csub{ilev}(Csub_rows,:);
          Caux = repmat ({Caux}, 1, rdim);
          Caux = blkdiag (Caux{:});

          dofs_on_lev = [];
          for icomp = 1:rdim
            dofs_on_lev = union (dofs_on_lev, (icomp-1)*hspace.ndof + (1:last_dof(ilev)));
          end
          A(dofs_on_lev,dofs_on_lev) = A(dofs_on_lev,dofs_on_lev) + Caux.' * A_side * Caux;

        end
      end
    end
  end

end
