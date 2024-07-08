% ADAPTIVITY_LINEAR_ELASTICITY_MIXED: solve the linear elasticity problem with an
% adaptive isogeometric method on a mixed formulation based on hierarchical Taylor-Hood splines.
%
% [geometry, hmsh, hspace, u, hspace_p, press, solution_data] = ...
%     adaptivity_linear_elasticity_mixed (problem_data, method_data, adaptivity_data, plot_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - symm_sides:   sides with symmetry boundary condition (may be empty)
%    - c_diff:       diffusion coefficient (see solve_laplace)
%    - grad_c_diff:  gradient of the diffusion coefficient (if not present, it is taken as zero)
%    - f:            function handle of the source term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
%
%  method_data : a structure with discretization data (see adaptivity_laplace). 
%     Using Taylor-Hood splines. The degree and regularity correspond to the Lagrange multiplier.
%  adaptivity_data: a structure with data for the adaptive method (see adaptivity_laplace)
%  plot_data: a structure to decide whether to plot things during refinement (see adaptivity_laplace)
%
% OUTPUT:
%    geometry:      geometry structure (see geo_load or mp_geo_load)
%    hmsh:          object representing the hierarchical mesh (see hierarchical_mesh)
%    hspace:        object representing the space of hierarchical splines for displacement (see hierarchical_space)
%    u:             computed degrees of freedom for the displacement, at the last iteration.
%    hspace_p:      object representing the space of hierarchical splines for the Lagrange multiplier (see hierarchical_space)
%    press:         computed degrees of freedom for the Lagrange multiplier, at the last iteration.
%    solution_data: output_data, see adaptivity_laplace.
% 
% Copyright (C) 2017-2018 Cesare Bracco
% Copyright (C) 2017-2024 Rafael Vazquez
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
% ADAPTIVITY_linear_elasticity: solve the linear elasticity problem with an adaptive isogeometric method based on hierarchical splines (mixed formulation).
%

function [geometry, hmsh, hspace, u, hspace_press, press, solution_data] = adaptivity_linear_elasticity_mixed (problem_data, method_data, adaptivity_data, plot_data)


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
nel = zeros (1, adaptivity_data.num_max_iter); ndof = nel; gest = nel+1; ndof2 = ndof;


% Initialization of the hierarchical mesh and spaces
[hmsh, hspace, hspace_press, geometry] = adaptivity_initialize_TH (problem_data, method_data);


% ADAPTIVE LOOP
iter = 0;
while (iter < adaptivity_data.num_max_iter)
  iter = iter + 1;
  
  if (plot_data.print_info)
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',iter);
  end

% SOLVE AND PLOT
  if (plot_data.print_info)
    disp('SOLVE:')
    fprintf('Number of elements: %d. Number of levels: %d \n', hmsh.nel, hspace.nlevels)
    fprintf('Number of DOFs (space 1): %d. Number of DOFs (space 2): %d \n', hspace.ndof, hspace_press.ndof);
  end
  [u,press] = adaptivity_solve_linear_elasticity_mixed (hmsh, hspace, hspace_press, problem_data);  %add outputs for condition number, sparsity e bandwidth?
  nel(iter) = hmsh.nel; ndof(iter) = hspace.ndof;  ndof2(iter) = hspace_press.ndof;
  
  if (plot_data.plot_hmesh)
    hmsh_plot_cells (hmsh, 10, fig_mesh);
    drawnow
  end
  if (plot_data.plot_discrete_sol)
    figure(fig_sol)  
    pts = {linspace(0, 1, 21), linspace(0, 1, 21)};  
    [eu, F] = sp_eval (u, hspace, geometry, pts);
    [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
    subplot (1,2,1)
    quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
    title ('Numerical solution'), axis equal tight
    subplot (1,2,2)
    eu2 = problem_data.uex (X, Y);
    quiver (X, Y, squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
    title ('Exact solution'), axis equal tight
    drawnow
  end
    
% ESTIMATE
  if (plot_data.print_info); disp('ESTIMATE:'); end
  est = adaptivity_estimate_linear_el_mixed(u, press, hmsh, hspace, hspace_press, problem_data, adaptivity_data);
  est_elem{iter} = est.';
  gest(iter) = norm (est);
  if (plot_data.print_info); fprintf('Computed estimate: %e \n', gest(iter)); end
  if (isfield (problem_data, 'graduex'))
    [~, err_l2(iter), err_h1s(iter),~,errl2_elem{iter},err_h1s_elem{iter}] = sp_h1_error (hspace, hmsh, u, problem_data.uex, problem_data.graduex);
    if (plot_data.print_info); fprintf('Error in H1 seminorm = %e\n', err_h1s(iter)); end
  elseif (isfield (problem_data, 'uex'))
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
%  marked{end} = union (marked{end}, 1);
  if (plot_data.print_info)
    fprintf('%d %s marked for refinement \n', num_marked, adaptivity_data.flag);
    disp('REFINE:')
  end

% REFINE
  [hmsh, hspace_aux] = adaptivity_refine_mixed (hmsh, {hspace, hspace_press}, marked, adaptivity_data);
  hspace = hspace_aux{1}; hspace_press = hspace_aux{2};
  clear hspace_aux
  fprintf('\n');
end

solution_data.iter = iter;
solution_data.gest = gest(1:iter);
solution_data.ndof = ndof(1:iter);
solution_data.ndof2 = ndof2(1:iter);
solution_data.nel  = nel(1:iter);

end
