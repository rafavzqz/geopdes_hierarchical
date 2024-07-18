% OP_PENALTY_DUDN: matrix and right-hand side to compute penalty terms to impose du/dn = f.
%  It computes the terms of the form Cp*(du/dn, dv/dn) and Cp*(f, dv/dn), 
%  using the same space for trial and test functions.
%
%   [mat, rhs] = op_penalty_dudn (hspace, hmsh, sides, Cpen, [coeff]);
%
% INPUT:
%
%  hspace: object representing the space of trial functions (see hierarchical_space_mp_C1)
%  hmsh:   object defining the domain partition and the quadrature rule (see hierarchical_mesh_mp)
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

function [A, rhs] = op_penalty_dudn (hspace, hmsh, nmnn_sides, Cpen, coeff)

  A = spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);
  rhs = zeros (hspace.ndof, 1);

  boundaries = hmsh.mesh_of_level(1).boundaries;
  Nbnd = cumsum ([0, boundaries.nsides]);
  last_dof = cumsum (hspace.ndof_per_level);

  for iref = nmnn_sides
    iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    for ilev = 1:hmsh.boundary.nlevels
      patch_bnd_shifting = cumsum ([0 hmsh.boundary.mesh_of_level(ilev).nel_per_patch]);
%      patch_shifting = cumsum ([0 hmsh.mesh_of_level(ilev).nel_per_patch]);
      dofs_on_lev = 1:last_dof(ilev);

      for ii = 1:numel(iref_patch_list)
        iptc_bnd = iref_patch_list(ii);
        iptc = boundaries(iref).patches(ii);
        iside = boundaries(iref).faces(ii);
        elems_patch = patch_bnd_shifting(iptc_bnd)+1:patch_bnd_shifting(iptc_bnd+1);
        [~, ~, elements] = intersect (hmsh.boundary.active{ilev}, elems_patch);

        if (~isempty (elements))
          msh_patch = hmsh.mesh_of_level(ilev).msh_patch{iptc};

          msh_side = msh_eval_boundary_side (msh_patch, iside, elements);
          msh_side_from_interior = msh_boundary_side_from_interior (msh_patch, iside);

          sp_bnd = hspace.space_of_level(ilev).sp_patch{iptc}.constructor (msh_side_from_interior);
          msh_side_from_interior_struct = msh_evaluate_element_list (msh_side_from_interior, elements);
          sp_bnd = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', true, 'gradient', true);

          % For simplicity I replaced charlen by 2^(-level)
%           coeff_at_qnodes =  coeff_at_qnodes * Cpen ./ msh_side.charlen;
          coeff_at_qnodes = Cpen * ones(msh_side.nqn, msh_side.nel)/ 2^(-ilev);

          [Cpatch, Cpatch_cols_lev] = sp_compute_Cpatch (hspace.space_of_level(ilev), iptc);
          [~,Csub_rows,Cpatch_cols] = intersect (hspace.Csub_row_indices{ilev}, Cpatch_cols_lev);
          Caux = Cpatch(:,Cpatch_cols) * hspace.Csub{ilev}(Csub_rows,:);

          tmp = op_gradu_n_gradv_n(sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
          tmp = Caux.' * tmp * Caux;

          A(dofs_on_lev,dofs_on_lev) = A(dofs_on_lev, dofs_on_lev) + tmp;

          if (nargin == 5)
            x = cell (hmsh.rdim, 1);
            for idim = 1:hmsh.rdim
              x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
            end
            rhs(dofs_on_lev) = rhs(dofs_on_lev) + ...
              Caux.' *op_gradv_n_f (sp_bnd, msh_side, Cpen*coeff(x{:}));
          end

        end
      end
    end
  end

end

