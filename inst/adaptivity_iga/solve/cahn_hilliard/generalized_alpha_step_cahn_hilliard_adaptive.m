% GENERALIZED_ALPHA_STEP_CAHN_HILLIARD_ADAPTIVE: perform one step of the generalized alpha method
%  for the solution of the Cahn-Hilliard equation, solving the nonlinear equation with Newton's method.
%  It is called from solve_step_adaptive_cahn_hilliard.
%
% INPUT:
%
%  u_n:          field at the previous time step.
%  u_dotn:       time derivative at the previous time step.
%  dt:           time step size
%  a_m:          parameter for the generalized alpha method
%  a_f:          parameter for the generalized alpha method
%  gamma:        parameter for the generalized alpha method
%  lambda:       parameter from the Cahn-Hilliard equation
%  mu, dmu:      function handle for value and derivative of the mu coefficient
%  Cpen:         penalization parameter for Nitsche's method
%  hspace:       space object (see hierarchical_space or hierarchical_space_mp_C1)
%  hmsh:         mesh object (see hierarchical_mesh or hierarchical_mesh_mp)
%  old_space:    matrices from the space of previous iterations. If not changed, some matrices are not recomputed.
%  nmnn_sides:   sides where to impose the Neumann condition.
%
% OUTPUT:
%
%  u_n1:      field at the new step.
%  u_dotn1:   time derivative at the new step.
%  old_space: matrices from the space of the previous iteration, updated
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

function [u_n1, udot_n1, old_space] = generalized_alpha_step_cahn_hilliard_adaptive (u_n, udot_n, dt, a_m, a_f, gamma, ...
                                       lambda, mu, dmu, Cpen, hspace, hmsh, old_space, nmnn_sides)

  % Convergence criteria
  n_max_iter = 20;
  tol_rel_res = 1e-10;
  tol_abs_res = 1e-10;

  % Predictor step
  u_n1 = u_n;
  udot_n1 = (gamma-1)/gamma * udot_n; 

  % Newton loop
  for iter = 0:n_max_iter

    % Field at alpha level
    udot_a = udot_n + a_m *(udot_n1-udot_n);
    u_a = u_n + a_f *(u_n1-u_n);

    % Compute the residual (internal)
    [Res_gl, stiff_mat, mass_mat, old_space] = ...
      Res_K_cahn_hilliard_adaptive(hspace, hmsh, lambda, ...
      Cpen, u_a, udot_a, mu, dmu, old_space, nmnn_sides);

    % Convergence check
    if (iter == 0)
      norm_res_0 = norm(Res_gl);
    end
    norm_res = norm(Res_gl);

    if (norm_res/norm_res_0 < tol_rel_res) % relative tolerance
      disp(strcat('iteration n°=',num2str(iter)))
      disp(strcat('norm (abs) residual=',num2str(norm_res)))
      break
    end
    if (norm_res<tol_abs_res) % absolute tolerance
      disp(strcat('iteration n°=',num2str(iter)))
      disp(strcat('norm absolute residual=',num2str(norm_res)))
      break
    end
    if (iter == n_max_iter)
      disp(strcat('Newton reached the maximum number of iterations'))
      disp(strcat('norm residual=',num2str(norm_res)))
    end
    
    % Compute the update, and update the solution
    A_gl = a_m * mass_mat + a_f * gamma *dt * stiff_mat ; 
    d_udot = - A_gl\Res_gl;

    udot_n1 = udot_n1 + d_udot;
    u_n1 = u_n1 + gamma * dt* d_udot;
  end

end
