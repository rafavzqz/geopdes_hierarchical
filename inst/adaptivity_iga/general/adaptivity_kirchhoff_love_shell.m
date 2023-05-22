% ADAPTIVITY_linear_elasticity: solve the linear elasticity problem with an
% adaptive isogeometric method based on hierarchical splines.
%
% [geometry, hmsh, hspace, u, solution_data] = adaptivity_linear_elasticity (problem_data, method_data, adaptivity_data, plot_data)
%
% XXXXXXXXXXXXX TEXT IS MISSING
% 
% Copyright (C) 2017 Cesare Bracco, Rafael Vazquez
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

function [geometry, hmsh, hspace, u, solution_data] = adaptivity_kirchhoff_love_shell (problem_data, method_data, adaptivity_data, plot_data)

if (nargin == 3)
  plot_data = struct ('print_info', true, 'plot_hmesh', false, 'plot_discrete_sol', false);
end
if (~isfield (plot_data, 'print_info'))
  plot_data.print_info = true;
end
if (~isfield (plot_data, 'plot_hmesh'))
  plot_data.plot_hmesh = false;
end
if (~isfield (plot_data, 'plot_discrete_sol'))
  plot_data.plot_discrete_sol = false;
end

% Initialization of some auxiliary variables
if (plot_data.plot_hmesh)
  fig_mesh = figure;
end
if (plot_data.plot_discrete_sol)
  fig_sol = figure;
end
nel = zeros (1, adaptivity_data.num_max_iter); ndof = nel; gest = nel+1;

% Initialization of the hierarchical mesh and space
[hmsh, hspace, geometry] = adaptivity_initialize_vector (problem_data, method_data);

% ADAPTIVE LOOP
iter = 0;
while (1)
  iter = iter + 1;

  if (plot_data.print_info)
    fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',iter);
  end
  
% SOLVE AND PLOT
  if (plot_data.print_info)
    disp('SOLVE:')
    fprintf('Number of elements: %d. Total DOFs: %d. Number of levels: %d \n', hmsh.nel, hspace.ndof, hspace.nlevels);
  end
  [u, energy(iter)] = adaptivity_solve_kirchhoff_love_shell (hmsh, hspace, problem_data);
  nel(iter) = hmsh.nel; ndof(iter) = hspace.ndof;
    
  if (plot_data.plot_hmesh)
    fig_mesh = hmsh_plot_cells (hmsh, 10, fig_mesh);
    drawnow
  end
  if (plot_data.plot_discrete_sol)
    npts = 51 * ones (1, hmsh.ndim);
    fig_sol = plot_numerical_and_exact_solution (u, hspace, geometry, npts, problem_data.uex, fig_sol); 
    drawnow
  end
  
% ESTIMATE
  if (plot_data.print_info); disp('ESTIMATE:'); end
  est = adaptivity_bubble_estimator_KL_shell (u, hmsh, hspace, problem_data, adaptivity_data);
  gest(iter) = norm (est);
  if (plot_data.print_info); fprintf('Computed error estimate: %e \n', gest(iter)); end
  if (isfield (problem_data, 'uex'))
    err_l2(iter) = sp_l2_error (hspace, hmsh, u, problem_data.uex);
    if (plot_data.print_info); fprintf('Error in L2 norm = %e\n', err_l2(iter)); end
  end
    
% STOPPING CRITERIA
  if (gest(iter) < adaptivity_data.tol)
    disp('Success: The error estimation reached the desired tolerance'); 
    solution_data.flag = 1; break
  elseif (iter == adaptivity_data.num_max_iter)
    disp('Warning: reached the maximum number of iterations')
    solution_data.flag = 2; break
  elseif (hmsh.nlevels >= adaptivity_data.max_level)
    disp('Warning: reached the maximum number of levels')
    solution_data.flag = 3; break
  elseif (hspace.ndof > adaptivity_data.max_ndof)
    disp('Warning: reached the maximum number of DOFs')
    solution_data.flag = 4; break
  elseif (hmsh.nel > adaptivity_data.max_nel)
    disp('Warning: reached the maximum number of elements')
    solution_data.flag = 5; break
  end
  
% MARK
  if (plot_data.print_info); disp('MARK:'); end
  [marked, num_marked] = adaptivity_mark (est, hmsh, hspace, adaptivity_data);
  if (plot_data.print_info)
    fprintf('%d %s marked for refinement \n', num_marked, adaptivity_data.flag);
    disp('REFINE:')
  end
  
% REFINE
  [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
  fprintf('\n');
end

solution_data.iter = iter;
solution_data.gest = gest(1:iter);
solution_data.ndof = ndof(1:iter);
solution_data.nel  = nel(1:iter);
if (exist ('err_l2', 'var'))
  solution_data.err_l2 = err_l2(1:iter);
end
solution_data.energy = energy(1:iter);

end