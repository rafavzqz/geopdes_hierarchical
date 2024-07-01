% ADAPTIVITY_KIRCHHOFF_LOVE_SHELL_MP_C1: solve the Kirchhoff-Love shell
%  problem in a multipatch domain with C^1 splines
%
% [geometry, hmsh, hspace, u, solution_data] = adaptivity_kirchhoff_love_shell_mp_C1 (problem_data, method_data, adaptivity_data, plot_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - drchlt_sides: sides with Dirichlet boundary condition (for now, only homogeneous conditions in every direction)
%    - E_coeff:      function handle for Young's modulus
%    - nu_coeff:     function handle for Poisson's ratio
%    - thickness:    scalar value, thickness of the shell
%    - f:            source term, distributed load
% 
%  method_data : a structure with discretization data (see adaptivity_laplace_mp_C1).
%  adaptivity_data: a structure with data for the adaptive method (see adaptivity_laplace_mp_C1).
%  plot_data: a structure to decide whether to plot things during refinement (see adaptivity_laplace_mp_C1).
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
%      - flag:       a flag with the stopping criterion satisfied (see adaptivity_laplace_mp_C1).
%
% The number of degrees of freedom is 3*hspace.ndof, that is, hspace only
%  contains scalar valued functions. This is different from the standard
%  implementation in GeoPDEs with the sp_vector (or sp_multipatch) class.
%
% Copyright (C) 2022-2024 Rafael Vazquez
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

function [geometry, hmsh, hspace, u, solution_data] = adaptivity_kirchhoff_love_shell_mp_C1 (problem_data, method_data, adaptivity_data, plot_data)

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
[hmsh, hspace, geometry] = adaptivity_initialize_laplace_mp_C1 (problem_data, method_data);

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
  [u, energy(iter)] = adaptivity_solve_kirchhoff_love_shell_mp_C1 (hmsh, hspace, problem_data, method_data);
  nel(iter) = hmsh.nel; ndof(iter) = hspace.ndof;

  if (plot_data.plot_hmesh)
    fig_mesh = hmsh_plot_cells (hmsh, 10, fig_mesh);
    drawnow
  end
  % if (plot_data.plot_discrete_sol)
  %   npts = 51 * ones (1, hmsh.ndim);
  %   fig_sol = plot_numerical_and_exact_solution (u, hspace, geometry, npts, problem_data.uex, fig_sol); 
  %   drawnow
  % end
  % if (isfield (plot_data, 'pts'))
  %   displ{iter} = sp_eval_phys (u, hspace, hmsh, geometry, plot_data.pts, plot_data.patch_for_pts);
  % end

% ESTIMATE
  if (plot_data.print_info); disp('ESTIMATE:'); end
  est = adaptivity_bubble_estimator_KL_shell (u, hmsh, hspace, problem_data, method_data, adaptivity_data);
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

if (isfield (plot_data, 'pts'))
  solution_data.displ = displ;
end

end