function [geometry, hmsh, hspace, u] = adaptivity_linear_elasticity (problem_data, method_data, adaptivity_data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the hierarchical mesh and space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hmsh, hspace, geometry] = adaptivity_initialize_vector (problem_data, method_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADAPTIVE LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter = 0;
while (iter < adaptivity_data.num_max_iter)
  iter = iter + 1;

  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',iter);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% SOLVE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('SOLVE:')
  u = adaptivity_solve_linear_elasticity (hmsh, hspace, problem_data);
  fprintf('Number of elements: %d. Total DOFs: %d \n', hmsh.nel, hspace.ndof);
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% ESTIMATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   disp('ESTIMATE:')
% %   est = estimate_laplace_param (u, hmsh, hspace, problem_data, adaptivity_data);
% %   gest = norm (est);
% %   fprintf('Computed estimate: %f \n', gest);
  if (isfield (problem_data, 'uex'))
    err_l2 = hspace_l2_error (hspace, hmsh, u, problem_data.uex);
    fprintf('Error in L2 norm = %g\n', err_l2);
  end
    
%   if (gest < adaptivity_data.tol)
%     disp('Success: The error estimation reached the desired tolerance')
%     break
%   end
%     
%   if (iter == adaptivity_data.num_max_iter)
%     disp('Warning: Maximum amount of iterations reached')
%     break
%   end
%     
%   if (hmsh.nlevels >= adaptivity_data.max_level)
%     disp('Warning: Maximum amount of levels reached')
%     break
%   end
%     
%   if (hspace.ndof > adaptivity_data.max_ndof)
%     disp('Warning: Maximum allowed DOFs achieved')
%     break
%   end
%     
%   if (hmsh.nel > adaptivity_data.max_nel)
%     disp('Warning: Maximum allowed amount of elements achieved')
%     break
%   end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% MARK
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   disp('MARK:')
%   [marked, num_marked] = adaptivity_mark (est, hmsh, hspace, adaptivity_data);
%   fprintf('%d %s marked for refinement \n', num_marked, adaptivity_data.flag);
  
  if (iter == 1)
    marked{1} = 1:hspace.ndof;
  elseif (iter == 2)
    marked{1} = [];
    marked{2} = 1:hspace.ndof;
  elseif (iter == 3)
    marked{1} = [];
    marked{2} = [];
    marked{3} = 1:hspace.ndof;
  else
    marked = cell (hspace.nlevels, 1);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% REFINE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('REFINE:')
  [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
  fprintf('\n');
end
