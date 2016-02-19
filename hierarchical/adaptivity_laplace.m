% ADAPTIVITY_LAPLACE: solve the Laplace problem with an adaptive isogeometric method based on hierarchical splines.
%
% [geometry, hmsh, hspace, u] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_diff:       diffusion coefficient (see solve_laplace)
%    - f:            source term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
%
%  method_data : a structure with discretization data. It contains the fields:
%    - degree:      degree of the spline functions.
%    - regularity:  continuity of the spline functions.
%    - nsub_coarse: number of subelements with respect to the geometry mesh (1 leaves the mesh unchanged)
%    - nsub_refine: number of subelements to be added at each refinement step (2 for dyadic)
%    - nquad:       number of points for Gaussian quadrature rule
%    - space_type:  'simplified' (only children of removed functions) or 'standard' (full hierarchical basis)
%    - truncated:   false (classical basis) or true (truncated basis)
%
%  adaptivity_data: a structure with data for the adaptive method. It contains the fields:
%    - flag:          refinement procedure, based either on 'elements' or on 'functions'
%    - mark_strategy: marking strategy. See 'adaptivity_mark' for details
%    - mark_param:    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX. See 'adaptivity_mark' for details
%    - max_level:     stopping criterium, maximum number of levels allowed during refinement
%    - max_ndof:      stopping criterium, maximum number of degrees of freedom allowed during refinement
%    - max_nel:       stopping criterium, maximum number of elements allowed during refinement
%    - num_max_iter:  stopping criterium, maximum number of iterations allowed
%    - tol:           stopping criterium, adaptive refinement is stopped when the norm of the estimator
%                      is lower than tol.
%
%  plot_data: a structure to decide whether to plot things during refinement.
%    - plot_hmesh
%    - plot_discrete_sol
%
% OUTPUT:
%    geometry: geometry structure (see geo_load)
%    hmsh:     object representing the hierarchical mesh (see hierarchical_mesh)
%    hspace:   object representing the space of hierarchical splines (see hierarchical_space)
%    u:        computed degrees of freedom, at the last iteration.
%    solution_data: a structure with the following fields
%      - iter:    iteration on which the adaptive procedure stopped
%      - ndof:    number of degrees of freedom for each computed iteration
%      - nel:     number of elements for each computed iteration
%      - gest:    norm of the estimator, for each computed iteration
%      - err_h1s: error in H1 seminorm for each iteration, if the exact solution is known
%      - flag:    a flag with one of the following values:
%          -1: the coefficients for the partition of unity were wrong. This is probably caused by a bug.
%           1: convergence is reached, the norm of the estimator is lower than the given tolerance
%           2: maximum number of iterations reached before convergence.
%           3: maximum number of levels reached before convergence
%           4: maximum number of degrees of freedom reached before convergence
%           5: maximum number of elements reached before convergence
% 
%    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX: anything else?
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
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


function  [geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data)

% Initialization of some auxiliary variables
if (plot_data.plot_hmesh)
  fig_mesh = figure;
end
if (plot_data.plot_discrete_sol)
  fig_sol = figure;
end
nel = zeros (1, adaptivity_data.num_max_iter); ndof = nel; gest = nel+1;

% Initialization of the hierarchical mesh and space
[hmsh, hspace, geometry] = adaptivity_initialize_laplace (problem_data, method_data);


% ADAPTIVE LOOP
iter = 0;
bpx = struct ('Ai', [], 'rhs', [], 'int_dofs', [], 'ndof', [], 'new_dofs', [], 'Pi', [], 'Qi', []);

bpx(1).Qi = speye (hspace.ndof);

while (iter < adaptivity_data.num_max_iter)
  iter = iter + 1;
  
  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',iter);
    
  if (~hspace_check_partition_of_unity (hspace, hmsh))
    disp('ERROR: The partition-of-the-unity property does not hold.')
    solution_data.flag = -1; break
  end

% SOLVE AND PLOT
  disp('SOLVE:')
  [u, bpx] = adaptivity_solve_laplace (hmsh, hspace, problem_data, bpx);
  fprintf('Number of elements: %d. Total DOFs: %d \n', hmsh.nel, hspace.ndof);
  nel(iter) = hmsh.nel; ndof(iter) = hspace.ndof;

  if (plot_data.plot_hmesh)
    hmsh_plot_cells (hmsh, 10, fig_mesh);
  end
  if (plot_data.plot_discrete_sol)
    figure(fig_sol)
    npts = 51 * ones (1, hmsh.ndim);
    plot_numerical_and_exact_solution (u, hspace, geometry, npts, problem_data.uex); 
  end
  if (plot_data.plot_hmesh || plot_data.plot_discrete_sol)
   disp('Paused. Type "dbcont" to continue')
   keyboard
  end

% ESTIMATE
  disp('ESTIMATE:')
  est = adaptivity_estimate_laplace (u, hmsh, hspace, problem_data, adaptivity_data);
  gest(iter) = norm (est);
  fprintf('Computed estimate: %f \n', gest(iter));
  if (isfield (problem_data, 'graduex'))
    [~, ~, err_h1s(iter)] = sp_h1_error (hspace, hmsh, u, problem_data.uex, problem_data.graduex);
    fprintf('Error in H1 seminorm = %g\n', err_h1s(iter));
  end

% STOPPING CRITERIA
  if (gest(iter) < adaptivity_data.tol)
    disp('Success: The error estimation reached the desired tolerance'); 
    solution_data.flag = 1; break
  elseif (iter == adaptivity_data.num_max_iter)
    disp('Warning: Maximum amount of iterations reached')
    solution_data.flag = 2; break
  elseif (hmsh.nlevels >= adaptivity_data.max_level)
    disp('Warning: Maximum amount of levels reached')
    solution_data.flag = 3; break
  elseif (hspace.ndof > adaptivity_data.max_ndof)
    disp('Warning: Maximum allowed DOFs achieved')
    solution_data.flag = 4; break
  elseif (hmsh.nel > adaptivity_data.max_nel)
    disp('Warning: Maximum allowed amount of elements achieved')
    solution_data.flag = 5; break
  end
  
% MARK
  disp('MARK:')
  [marked, num_marked] = adaptivity_mark (est, hmsh, hspace, adaptivity_data);
  fprintf('%d %s marked for refinement \n', num_marked, adaptivity_data.flag);

% REFINE
  disp('REFINE:')
  if (method_data.bpx_dofs)
    [hmsh, hspace, Cref] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
  else
    [hmsh, hspace, Cref, new_dofs] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
    % NEW_DOFS VA POI INTERSECATO CON int_dofs dentro solve
    bpx(iter+1).new_dofs = new_dofs;
  end
  fprintf('\n');
  
%
bpx(iter).Pi = Cref;
bpx(iter+1).Qi = speye (hspace.ndof);
for lind = iter:-1:1
  bpx(lind).Qi = bpx(lind+1).Qi * bpx(lind).Pi;
end

end


solution_data.iter = iter;
solution_data.gest = gest(1:iter);
solution_data.ndof = ndof(1:iter);
solution_data.nel  = nel(1:iter);
if (exist ('err_h1s', 'var'))
  solution_data.err_h1s = err_h1s(1:iter);
end

end