% ADAPTIVITY_BILAPLACE_MP_C1: solve the bilaplacian problem with an adaptive 
%  isogeometric method based on hierarchical splines for C^1 multipatch splines.
%
% [geometry, hmsh, hspace, u, solution_data] = adaptivity_bilaplace (problem_data, method_data, adaptivity_data, plot_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - drchlt_sides: sides with essential boundary conditions (value and normal derivative)
%    - c_diff:       physical parameter (see adaptivity_solve_bilaplace_mp_C1)
%    - f:            function handle of the source term
%    - g:            function for Neumann condition
%    - h:            function for Dirichlet boundary condition
%    - uex, graduex: exact solution. If known, the boundary conditions are computed from it.
%
%  method_data : a structure with discretization data. It contains the fields:
%    - degree:       degree of the spline functions.
%    - regularity:   continuity of the spline functions.
%    - nsub_coarse:  number of subelements with respect to the geometry mesh (1 leaves the mesh unchanged)
%    - nsub_refine:  number of subelements to be added at each refinement step (2 for dyadic)
%    - nquad:        number of points for Gaussian quadrature rule
%    - space_type:   'simplified' (only children of removed functions) or 'standard' (full hierarchical basis)
%    - truncated:    false (classical basis) or true (truncated basis)
%
%  adaptivity_data: a structure with data for the adaptive method. It contains the fields:
%    - flag:          refinement procedure, based either on 'elements' or on 'functions'
%    - mark_strategy: marking strategy. See 'adaptivity_mark' for details
%    - mark_param:    a parameter to decide how many entities should be marked. See 'adaptivity_mark' for details
%    - max_level:     stopping criterium, maximum number of levels allowed during refinement
%    - max_ndof:      stopping criterium, maximum number of degrees of freedom allowed during refinement
%    - max_nel:       stopping criterium, maximum number of elements allowed during refinement
%    - num_max_iter:  stopping criterium, maximum number of iterations allowed
%    - tol:           stopping criterium, adaptive refinement is stopped when the global error estimator
%                      is lower than tol.
%    - C0_est:        an optional multiplicative constant for scaling the error estimators (default value: 1).
%    - adm_class:     admissibility class, to control the interaction of functions of different levels;
%    - adm_type:      either 'T-admissible' or 'H-admissible'
%
%  plot_data: a structure to decide whether to plot things during refinement.
%    - plot_hmesh:        plot the mesh at every iteration
%    - plot_discrete_sol: plot the discrete solution at every iteration
%    - print_info:        display info on the screen on every iteration (number of elements, 
%                          number of functions, estimated error, number of marked elements/functions...)
%
% OUTPUT:
%    geometry:      geometry structure (see mp_geo_load)
%    hmsh:          object representing the hierarchical mesh (see hierarchical_mesh_mp)
%    hspace:        object representing the space of hierarchical splines (see hierarchical_space_mp_C1)
%    u:             computed degrees of freedom, at the last iteration.
%    solution_data: a structure with the following fields
%      - iter:       iteration on which the adaptive procedure stopped
%      - ndof:       number of degrees of freedom for each computed iteration
%      - nel:        number of elements for each computed iteration
%      - gest:       global error estimator, for each computed iteration
%      - err_l2:     error in L2 norm for each iteration, if the exact solution is known
%      - err_h1s:    error in H1 seminorm for each iteration, if the exact solution is known
%      - err_h1:     error in H1 norm for each iteration, if the exact solution is known
%      - err_h2s:    error in H2 seminorm (only Laplacian for surfaces) for each iteration, if the exact solution is known
%      - err_h2:     error in H2 norm (only Laplacian for surfaces) for each iteration, if the exact solution is known
%      - flag:       a flag with one of the following values:
%          -1: the coefficients for the partition of unity were wrong. This is probably caused by a bug.
%           1: convergence is reached, the global estimator is lower than the given tolerance
%           2: maximum number of iterations reached before convergence.
%           3: maximum number of levels reached before convergence
%           4: maximum number of degrees of freedom reached before convergence
%           5: maximum number of elements reached before convergence
%
%
% Copyright (C) 2022-2023 Rafael Vazquez
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


function  [geometry, hmsh, hspace, u, solution_data] = adaptivity_bilaplace_mp_C1 (problem_data, method_data, adaptivity_data, plot_data)

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
if (~isfield (problem_data, 'uex'))
  problem_data.uex = [];
end

% Initialization of some auxiliary variables
if (plot_data.plot_hmesh)
  fig_mesh = figure;
end
if (plot_data.plot_discrete_sol)
  fig_sol = figure;
end
nel = zeros (1, adaptivity_data.num_max_iter); ndof = nel; gest = nel+NaN;

if (isfield (problem_data, 'hessuex'))
    err_h2 = gest;
    err_h1 = gest;
    err_l2 = gest;
    err_h2s = gest;
    err_h1s = gest;
elseif (isfield (problem_data, 'graduex'))
    err_h1 = gest;
    err_l2 = gest;
    err_h1s = gest;
end
  
% Initialization of the hierarchical mesh and space
if (~isfield(method_data, 'interface_regularity') || method_data.interface_regularity ~= 1)
  warning('Setting interface regularity to C1')
  method_data.interface_regularity = 1;
end
[hmsh, hspace, geometry] = adaptivity_initialize_laplace (problem_data, method_data);

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
  u = adaptivity_solve_bilaplace_mp_C1 (hmsh, hspace, problem_data);
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
  est = adaptivity_bubble_estimator_bilaplace (u, hmsh, hspace, problem_data, adaptivity_data);
%  est = adaptivity_estimate_bilaplace_h_h2 (u, hmsh, hspace, problem_data);
  gest(iter) = norm (est);
  if (plot_data.print_info); fprintf('Computed error estimate: %f \n', gest(iter)); end
  if (isfield (problem_data, 'hessuex'))
    [err_h2(iter), err_h1(iter), err_l2(iter), err_h2s(iter), err_h1s(iter)] = sp_h2_error (hspace, hmsh, u, problem_data.uex, problem_data.graduex, problem_data.hessuex);
    if (plot_data.print_info); fprintf('Error in H2 seminorm = %g\n', err_h2s(iter)); end
  elseif (isfield (problem_data, 'lapuex'))
    [err_h2(iter), err_h1(iter), err_l2(iter), err_h2s(iter), err_h1s(iter)] = sp_h2_equiv_lap_error (hspace, hmsh, u, problem_data.uex, problem_data.graduex, problem_data.lapuex);
    if (plot_data.print_info); fprintf('Error in (equivalent) H2 seminorm = %g\n', err_h2s(iter)); end
  elseif (isfield (problem_data, 'graduex'))
    [err_h1(iter), err_l2(iter), err_h1s(iter)] = sp_h1_error (hspace, hmsh, u, problem_data.uex, problem_data.graduex);
    if (plot_data.print_info); fprintf('Error in H1 seminorm = %g\n', err_h1s(iter)); end
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
end

solution_data.iter = iter;
solution_data.gest = gest(1:iter);
solution_data.ndof = ndof(1:iter);
solution_data.nel  = nel(1:iter);
if (exist ('err_h2s', 'var'))
  solution_data.err_h2s = err_h2s(1:iter);
  solution_data.err_h2 = err_h2(1:iter);
end
if (exist ('err_h1s', 'var'))
  solution_data.err_h1s = err_h1s(1:iter);
  solution_data.err_h1 = err_h1(1:iter);
  solution_data.err_l2 = err_l2(1:iter);
end

end