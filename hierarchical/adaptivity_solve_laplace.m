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

if (adaptivity_data.num_max_iter <= 0)
  return
end

while 1
    iter = iter + 1;
    
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n',iter);
    
    if (~check_partition_of_the_unity(hmsh, hspace))
        disp('ERROR: The partition-of-the-unity property does not hold.')
        return,
    end
    
    if (plot_hmesh)
      hmsh_plot_cells (hmsh, 1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SOLVE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    u = adaptivity_solve (hmsh, hspace, problem_data);
    
    if (plot_discrete_sol)
       plot_numerical_and_exact_solution (u, hmsh, hspace, problem_data.uex); 
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Error computation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [err_h1, err, err_h1s] = hsp_h1_error (hspace, hmsh, u, problem_data.uex, problem_data.graduex);
    fprintf('error_H1s = %g\n', err_h1s);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ESTIMATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    disp('Computing error estimators:')
    est = estimate_laplace_param (u, hmsh, hspace, problem_data, adaptivity_data);
    gest = norm (est);
    tempo = toc;
    fprintf('EST: %f (%f seconds)\n', gest, tempo);
    
    if (gest < adaptivity_data.tol)
        disp('Success: The error estimation reached the desired tolerance')
        return,
    end
    
    if (iter == adaptivity_data.num_max_iter)
        disp('Warning: Maximum amount of iterations reached')
        return;
    end
    
    if (hmsh.nlevels >= adaptivity_data.max_level)
        disp('Warning: Maximum amount of levels reached')
        return;
    end
    
    if (hspace.ndof > adaptivity_data.max_ndof)
        disp('Warning: Maximum allowed DOFs achieved')
        return;
    end
    
    if (hmsh.nel > adaptivity_data.max_nel)
        disp('Warning: Maximum allowed amount of elements achieved')
        return;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MARK
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tic
    disp('Marking:')
    [marked, num_marked] = mark (est, hmsh, hspace, adaptivity_data);
    tempo = toc;
    fprintf('%d %s marked for refinement (%f seconds)\n', num_marked, adaptivity_data.flag, tempo);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% REFINE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [hmsh, hspace] = refine (hmsh, hspace, marked, adaptivity_data.flag);
    
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');
end
