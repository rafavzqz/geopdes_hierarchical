function  [geometry, hmsh, hspace, u, gest, err_h1s, iter] = adaptivity_solve_laplace(problem_data, method_data, adaptivity_data, plot_hmesh, plot_discrete_sol)
%
% function  [hmsh, hspace, u, gest, err_h1s, iter] = adaptivity_solve_laplace(problem_data, method_data, adaptivity_data, plot_hmesh, plot_discrete_sol)
%
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% ATENCION: Mejorar esta funcion, decidir bien los argumentos de entrada y
% salida
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the hierarchical mesh and space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hmsh, hspace, geometry] = adaptivity_initialize (problem_data, method_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADAPTIVE LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter = 0;
while (iter < adaptivity_data.num_max_iter)
  iter = iter + 1;
  
  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',iter);
    
  if (~hspace_check_partition_of_unity (hspace, hmsh))
    disp('ERROR: The partition-of-the-unity property does not hold.')
    break
  end
    
  if (plot_hmesh)
    hmsh_plot_cells (hmsh, 1);
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% SOLVE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('SOLVE:')
  u = adaptivity_solve (hmsh, hspace, problem_data);
  fprintf('Number of elements: %d. Total DOFs: %d \n', hmsh.nel, hspace.ndof);
    
  if (plot_discrete_sol)
    plot_numerical_and_exact_solution (u, hmsh, hspace, problem_data.uex); 
  end
       
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% ESTIMATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('ESTIMATE:')
  est = estimate_laplace_param (u, hmsh, hspace, problem_data, adaptivity_data);
  gest = norm (est);
  fprintf('Computed estimate: %f \n', gest);
  if (isfield (problem_data, 'graduex'))
    [err_h1, err_l2, err_h1s] = hspace_h1_error (hspace, hmsh, u, problem_data.uex, problem_data.graduex);
    fprintf('Error in H1 seminorm = %g\n', err_h1s);
  end
    
  if (gest < adaptivity_data.tol)
    disp('Success: The error estimation reached the desired tolerance')
    break
  end
    
  if (iter == adaptivity_data.num_max_iter)
    disp('Warning: Maximum amount of iterations reached')
    break
  end
    
  if (hmsh.nlevels >= adaptivity_data.max_level)
    disp('Warning: Maximum amount of levels reached')
    break
  end
    
  if (hspace.ndof > adaptivity_data.max_ndof)
    disp('Warning: Maximum allowed DOFs achieved')
    break
  end
    
  if (hmsh.nel > adaptivity_data.max_nel)
    disp('Warning: Maximum allowed amount of elements achieved')
    break
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% MARK
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('MARK:')
  [marked, num_marked] = adaptivity_mark (est, hmsh, hspace, adaptivity_data);
  fprintf('%d %s marked for refinement \n', num_marked, adaptivity_data.flag);
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% REFINE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('REFINE:')
  [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
  fprintf('\n');
end
