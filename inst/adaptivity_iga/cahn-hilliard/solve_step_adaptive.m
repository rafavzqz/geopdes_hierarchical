%--------------------------------------------------------------------------
% adaptivity in space
%--------------------------------------------------------------------------
function [u_n1, udot_n1, hspace, hmsh, est, old_space] = solve_step_adaptive...
  (u_n, udot_n, hspace, hmsh, dt, a_m, a_f, gamma, Cpen, ...
   problem_data, adaptivity_data, old_space, nmnn_sides)

  lambda = problem_data.lambda;
  mu = problem_data.mu;
  dmu = problem_data.dmu;
  iter = 0;
  while (1)
    iter = iter + 1;
    disp(strcat('%%%%%%%%%%%%%%%%% Adaptivity iteration ',num2str(iter),' %%%%%%%%%%%%%%%%%'));
        
    %------------------------------------------------------------------
    % solve
    [u_n1, udot_n1, old_space] = generalized_alpha_step(u_n, udot_n, dt, a_m, a_f, gamma, lambda, mu, dmu, ...
                                                        Cpen, hspace, hmsh, old_space, nmnn_sides);

    %------------------------------------------------------------------
    %estimate
% %%%%%%%%%%%%%%%%%%%%%%%%%% TODO: write the help of the function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    est = adaptivity_estimate_cahn_hilliard (u_n1, hmsh, hspace, adaptivity_data);

    %------------------------------------------------------------------
    % stopping criteria
    if (iter == adaptivity_data.num_max_iter)
      disp('Warning: reached the maximum number of iterations')
    elseif (hmsh.nlevels > adaptivity_data.max_level)
      disp(strcat('number of levels =',num2str(hmsh.nlevels))) 
      disp('Warning: reached the maximum number of levels')
      break
    end
        
    if (hmsh.nlevels == adaptivity_data.max_level && hmsh.nel == hmsh.nel_per_level(end))
      disp('Mesh completely refined at the maximum level')
      break
    end

    %------------------------------------------------------------------
    % mark
% %%%%%%%%%%%%%%%%%%%%%%%%%% TODO: write the help of the function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    [marked, num_marked_ref] = adaptivity_mark_cahn_hilliard (est, hmsh, hspace, adaptivity_data);
    % limit the maximum refinement depth
    if (hmsh.nlevels == adaptivity_data.max_level)
      num_deleted = numel(marked{hmsh.nlevels});
      marked{hmsh.nlevels} = [];
      num_marked_ref = num_marked_ref - num_deleted ;
    end

    % stopping criterion
    if (num_marked_ref == 0)
      disp('No element is refined')
      break  
    end

    % refine
    [hmsh, hspace, Cref] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
    old_space = struct ('modified', true, 'space', [], 'mesh', [], 'mass_mat', [], ...
                        'lapl_mat', [], 'bnd_mat', [], 'Pen', [], 'pen_rhs', []);

    % recompute control variables
    u_n = Cref * u_n;       
    udot_n = Cref * udot_n;
   
  end % end loop adaptivity
end
