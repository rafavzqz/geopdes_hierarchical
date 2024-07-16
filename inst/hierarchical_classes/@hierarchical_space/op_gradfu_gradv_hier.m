% OP_GRADFU_GRADV_HIER: assemble the matrix A = [a(i,j)], a(i,j) = (grad (f(w) u_j), grad v_i) = 
%  = (f(w) grad u_j, grad v_i) + (u_j grad (f(w)), grad v_i), 
%  for hierarchical splines, using the multilevel structure, with "f" a given function
%  (derivative also needed), and "w" a discrete solution defined on the
%  same space. Useful to compute nonlinear terms.
%
%   [mat1, mat2] = op_gradfu_gradv_hier (hspace, hmsh, uhat, f, df);
%
% INPUT:
%
%   hspace:  object representing the hierarchical space of test functions  (see hierarchical_space)
%   hmsh:    object representing the hierarchical mesh (see hierarchical_mesh)
%   uhat:    degrees of freedom of "w", of size hspace.ndof x 1.
%   f:       function handle to compute f(w).
%   df:      function handle to compute df(w), the gradient of f.
%
% OUTPUT:
%
%   mat1:   assembled matrix for (f(w) grad u_j, grad v_i)
%   mat2:   assembled matrix for (u_j grad (f(w)), grad v_i)
%
% Copyright (C) 2023-2024, Michele Torre, Rafael Vazquez
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

function [A, B] = op_gradfu_gradv_hier (hspace, hmsh, uhat, f, df)

  A = spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);
  B = spalloc (hspace.ndof, hspace.ndof, 6*hspace.ndof);

  ndofs_u = 0;
  ndofs_v = 0;
  for ilev = 1:hmsh.nlevels
    ndofs_u = ndofs_u + hspace.ndof_per_level(ilev);
    ndofs_v = ndofs_v + hspace.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)

      % space
      spu_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true, 'gradient', true, 'laplacian', true);
      spu_lev = change_connectivity_localized_Csub (spu_lev, hspace, ilev);

      % coefficients
      ndof_until_lev = sum (hspace.ndof_per_level(1:ilev));
      uhat_lev = hspace.Csub{ilev} * uhat(1:ndof_until_lev);

      utemp = sp_eval_msh (uhat_lev, spu_lev, hmsh.msh_lev{ilev}, {'value', 'gradient'});
      u = utemp{1};
      gradu = utemp{2};

      coeffs_A = f(u);
      coeffs_B = df(u);
      coeffs_Bv = gradu;
      for idim = 1:hmsh.ndim
        coeffs_Bv(idim,:,:)=coeffs_Bv(idim,:,:) .* reshape(coeffs_B, 1, size(coeffs_B,1), size(coeffs_B,2) );
      end

      % matrices
      A_lev = op_gradu_gradv (spu_lev, spu_lev, hmsh.msh_lev{ilev}, coeffs_A);
      B_lev = op_vel_dot_gradu_v (spu_lev, spu_lev, hmsh.msh_lev{ilev}, coeffs_Bv)';
        
      % assemble
      dofs_u = 1:ndofs_u;
      dofs_v = 1:ndofs_v;

      A(dofs_v,dofs_u) = A(dofs_v,dofs_u) + hspace.Csub{ilev}.' * A_lev * hspace.Csub{ilev};   
      B(dofs_v,dofs_u) = B(dofs_v,dofs_u) + hspace.Csub{ilev}.' * B_lev * hspace.Csub{ilev};   
    end
  end

end
