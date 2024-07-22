% RES_K_CAHN_HILLIARD: Cahn-Hilliard equation, computation of the residual
% for Newton's method, and of some necessary matrices for the method. 
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

function [Res_gl, stiff_mat, mass_mat, old_space] = Res_K_cahn_hilliard(hspace, hmsh, lambda, Cpen, u_a, udot_a, ...
                                                      mu, dmu, old_space, nmnn_sides)

  % Double well (matrices)
  [term2, term2K] = op_gradfu_gradv_hier (hspace, hmsh, u_a, mu, dmu);

  if (old_space.modified == true)
    
    % mass matrix
    mass_mat = op_u_v_hier (hspace,hspace,hmsh);

    % bilaplacian (matrix)
    lapl_mat = op_laplaceu_laplacev_hier (hspace, hspace, hmsh, lambda);

    % Compute the boundary term (Nitsche method). The right-hand side is set to zero.
    bnd_mat = op_nitsche_consistency_cahn_hilliard (hspace, hmsh, nmnn_sides, lambda);
    [Pen, pen_rhs] = op_penalty_dudn (hspace, hmsh, nmnn_sides, Cpen);

    % update old_space
    old_space = struct ('modified', false, ...
      'mass_mat', mass_mat, 'lapl_mat', lapl_mat, 'bnd_mat', bnd_mat, 'Pen', Pen, 'pen_rhs', pen_rhs);
      
  elseif (old_space.modified == false)
    mass_mat =  old_space.mass_mat;
    lapl_mat =  old_space.lapl_mat;
    bnd_mat =  old_space.bnd_mat;
    Pen = old_space.Pen;
    pen_rhs =  old_space.pen_rhs;
  end

  %----------------------------------------------------------------------
  % Residual
  Res_gl = mass_mat*udot_a + term2*u_a + lapl_mat*u_a - (bnd_mat + bnd_mat.')*u_a + Pen*u_a - pen_rhs;

  % Tangent stiffness matrix (mass is not considered here)
  stiff_mat = term2 + term2K + lapl_mat - (bnd_mat + bnd_mat.') + Pen;

end
