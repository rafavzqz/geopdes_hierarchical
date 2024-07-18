%--------------------------------------------------------------------------
% One step of generalized alpha method
%--------------------------------------------------------------------------
function [u_n1, udot_n1, old_space] = generalized_alpha_step(u_n, udot_n, dt, a_m, a_f, gamma, ...
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
      Res_K_cahn_hilliard(hspace, hmsh, lambda, Cpen, u_a, udot_a, mu, dmu, old_space, nmnn_sides);

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
