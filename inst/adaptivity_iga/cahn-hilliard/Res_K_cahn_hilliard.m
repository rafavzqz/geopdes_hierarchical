%--------------------------------------------------------------------------
% Cahn-Hilliard equation residual and tangent matrix
%--------------------------------------------------------------------------
function [Res_gl, stiff_mat, mass_mat, old_space] =  Res_K_cahn_hilliard(hspace, hmsh, lambda, Cpen, u_a, udot_a, ...
                                                      mu, dmu, old_space, nmnn_sides)

  if (old_space.modified == true)
    
    % mass matrix
    mass_mat = op_u_v_hier (hspace,hspace,hmsh);

    % Double well (matrices)
    [term2, term2K] = op_gradfu_gradv_hier (hspace, hmsh, u_a, mu, dmu);   

    % bilaplacian (matrix)
    lapl_mat = op_laplaceu_laplacev_hier (hspace, hspace, hmsh, lambda);

    % Compute the boundary term (Nitsche method). The right-hand side is set to zero.
    bnd_mat = int_boundary_term (hspace, hmsh, lambda, nmnn_sides);
    [Pen, pen_rhs] = op_penalty_dudn (hspace, hmsh, nmnn_sides, Cpen);

    % update old_space
    old_space = struct ('modified', false, 'space', hspace, 'mesh', hmsh, ...
      'mass_mat', mass_mat, 'lapl_mat', lapl_mat, 'bnd_mat', bnd_mat, 'Pen', Pen, 'pen_rhs', pen_rhs);
      
  elseif (old_space.modified == false)
    mass_mat =  old_space.mass_mat;
    lapl_mat =  old_space.lapl_mat;
    bnd_mat =  old_space.bnd_mat;
    Pen = old_space.Pen;
    pen_rhs =  old_space.pen_rhs;

    [term2, term2K] = op_gradfu_gradv_hier (old_space.space, old_space.mesh, u_a, mu, dmu);
  end

  %----------------------------------------------------------------------
  % Residual
  Res_gl = mass_mat*udot_a + term2*u_a + lapl_mat*u_a - (bnd_mat + bnd_mat.')*u_a + Pen*u_a - pen_rhs;

  % Tangent stiffness matrix (mass is not considered here)
  stiff_mat = term2 + term2K + lapl_mat - (bnd_mat + bnd_mat.') + Pen;

end
