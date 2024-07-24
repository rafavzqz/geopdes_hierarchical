% ADAPTIVITY_CAHN_HILLIARD_MP_C1: solve the Cahn-Hilliard equation, with a generalized alpha discretization in time, 
%  and adaptive (refining and coarsening) C1-multipatch hierarchical splines in space.
%
% The function solves the problem of finding u such that
%
%  du/dt - Delta (mu(u) - lambda*Delta u) = 0
%
% with Delta the Laplacian, and mu(u) = alpha u^3 - beta u, and Neumann boundary conditions.
%
% For details on the problem and the formulation, see
%  H. Gomez, V.M. Calo, Y. Bazilevs, T.J.R. Hughes, CMAME 197 (2008), 4333-4352.
%  H. Gomez, A. Reali, G. Sangalli, J. Comput. Physics 262 (2014), 153-171.
%  C. Bracco, C. Giannelli, A. Reali, M. Torre, R. Vazquez, CMAME 417 (2023), 116355. 
%
% USAGE:
%
%   [geometry, hmsh, hspace, results] = adaptivity_cahn_hilliard (problem_data, ...
%                  method_data, adaptivity_data, save_info)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - lambda:       parameter representing the length scale of the problem, and the width of the interface
%    - mu:           function handle to compute mu (from the double well function)
%    - dmu:          function handle to compute the derivative of mu
%    - initial_time: initial time of the simulation
%    - Time_max:     final time
%    - fun_u:        initial condition. Equal to zero by default.
%    - fun_udot:     initial condition for time derivative. Equal to zero by default.
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:      degree of the spline functions.
%    - regularity:  continuity of the spline functions.
%    - nsub_coarse: number of subelements with respect to the geometry mesh (1 leaves the mesh unchanged)
%    - nsub_refine: number of subelements to be added at each refinement step (2 for dyadic)
%    - nquad:       number of points for Gaussian quadrature rule
%    - space_type:  'simplified' (only children of removed functions) or 'standard' (full hierarchical basis)
%    - truncated:   false (classical basis) or true (truncated basis)
%    - interface_regularity: must be set to 1.
%    - dt:          time step size for generalized-alpha method
%    - rho_inf_gen_alpha: parameter in [0,1], which governs numerical damping of the generalized alpha method
%    - Cpen_projection: penalty parameter to impose zero flux at the initial condition
%    - Cpen_Nitsche:    penalty parameter for Nitsche's method
%
%  adaptivity_data: a structure with data for the adaptive method. It contains the fields:
%    - flag:          refinement procedure, based either on 'elements' or on 'functions'
%    - mark_strategy: marking strategy. See 'adaptivity_mark' for details
%    - mark_param:    a parameter to decide how many entities should be marked. See 'adaptivity_mark' for details
%    - max_level:     stopping criterium, maximum number of levels allowed during refinement
%    - num_max_iter:  stopping criterium, maximum number of iterations allowed
%    - estimator_type: either 'field' or 'gradient'
%    - adm_class:     admissibility class, to control the interaction of functions of different levels;
%    - adm_type:      either 'T-admissible' or 'H-admissible'
%    - time_delay:    time delay to activate coarsening
%
%  save_info: a structure with information about when and where to save the results
%    - folder_name: folder in which to save the results.
%    - file_name:   name of the files, a number will be appended to this.
%    - time_save:   time steps at which the solution should be saved.
%    - vtk_pts:     cell-array with the univariate points to export the VTK file (see sp_to_vtk)
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  hmsh:     hierarchical mesh object at the last computed step (see hierarchical_mesh)
%  hspace:   hierarchical space object at the last computed step (see hierarchical_space)
%  results:  a struct with the saved results, containing the following fields:
%    - time: (array of length Ntime) time at which the solution was saved
%    Since the mesh and space are changed during the simulation, to reduce the memory usage
%     they are saved in files corresponding to the time in results.time (see save_info)
%
% Only periodic and Neumann boundary conditions are implemented. Neumann
%  conditions are considered by default.
%
% Copyright (C) 2023, 2024 Michele Torre, Rafael Vazquez
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

function [geometry, hmsh, hspace, results] = ...
             adaptivity_cahn_hilliard_mp_C1 (problem_data, method_data, adaptivity_data, save_info)

%%-------------------------------------------------------------------------
% Initialization of the coarsest level of the hierarchical mesh and space
if (~isfield(method_data, 'interface_regularity') || method_data.interface_regularity ~= 1)
  warning('Setting interface regularity to C1')
  method_data.interface_regularity = 1;
end
[hmsh, hspace, geometry] = adaptivity_initialize_laplace (problem_data, method_data);

if (exist ('nmnn_sides','var') && ~isempty (nmnn_sides))
  disp('User defined Neumann sides deleted. Neumann conditions used everywhere.')
  clear nmnn_sides
end
nmnn_sides = 1:numel(hmsh.mesh_of_level(1).boundaries);

% Refine the mesh up to a predefined level
n_refinements = adaptivity_data.max_level - 1; % number of uniform refinements
[hmsh, hspace] = uniform_refinement(hmsh, hspace, n_refinements, adaptivity_data);

%%-------------------------------------------------------------------------
% Initial conditions, with a penalty term
mass_mat = op_u_v_hier (hspace,hspace,hmsh);
[Pen, ~] = op_penalty_dudn (hspace, hmsh, nmnn_sides, method_data.Cpen_projection);
mass_proj = mass_mat + Pen;

if (isfield(problem_data,'fun_u'))
  if (isnumeric(problem_data.fun_u))
    u_n = fun_u;
  else
    rhs = op_f_v_hier(hspace, hmsh, problem_data.fun_u);
    u_n = mass_proj \ rhs;
  end
else
  u_n = zeros(hspace.ndof,1);
end
    
if (isfield(problem_data,'fun_udot'))
  if (isnumeric(problem_data.fun_udot))
    udot_n = fun_udot;
  else
    rhs = op_f_v_hier(hspace, hmsh, problem_data.fun_udot);
    udot_n = mass_proj \ rhs;
  end
else
  udot_n = zeros(hspace.ndof,1);
end
clear mass_proj

%%-------------------------------------------------------------------------
% Generalized alpha parameters
rho_inf_gen_alpha = method_data.rho_inf_gen_alpha;
a_m = .5*((3-rho_inf_gen_alpha)/(1+rho_inf_gen_alpha));
a_f =  1/(1+rho_inf_gen_alpha);
gamma =  .5 + a_m - a_f;

%%-------------------------------------------------------------------------
% save matrices previous mesh
old_space = struct ('modified', true, 'mass_mat', [], ...
  'lapl_mat', [], 'bnd_mat', [], 'Pen', [], 'pen_rhs', []);

%%-------------------------------------------------------------------------
% Initialize structure to store the results
time = problem_data.initial_time;
time_save = save_info.time_save;
time_save = time_save(time_save>=problem_data.initial_time & time_save<=problem_data.Time_max);

results.time = zeros(length(time_save), 1);
time_save(end+1) = problem_data.Time_max + 1e5;

% Save initial conditions
save_id = 1;
if (time >= time_save(1))
  save_results_step(u_n, udot_n, time, hspace, hmsh, geometry, save_info, save_id)
  results.time(save_id) = time;
  save_id = 2;
  time_save = time_save(time_save > time);
end

%%-------------------------------------------------------------------------
% Loop over time steps
dt = method_data.dt;
while time < problem_data.Time_max

  disp('----------------------------------------------------------------')
  disp(strcat('time step t=',num2str(time)))
  disp(strcat('Number of elements = ', num2str(hmsh.nel)))

  %----------------------------------------------------------------------
  % one time step with adaptivity (refinement) in space
  [u_n1, udot_n1, hspace, hmsh, est, mark_param_coarsening, old_space] = solve_step_adaptive_cahn_hilliard ...
                              (u_n, udot_n, hspace, hmsh, dt, a_m, a_f, gamma, ...
                               method_data.Cpen_nitsche, problem_data, ...
                               adaptivity_data, old_space, nmnn_sides);
    
  %----------------------------------------------------------------------
  % coarsening
  adaptivity_data.mark_param_coarsening =  mark_param_coarsening;
  if (time >= adaptivity_data.time_delay)
    [u_n1, udot_n1, hspace, hmsh, old_space] = coarsening_algorithm ...
      (est, hmsh, hspace, adaptivity_data, u_n1, udot_n1, method_data.Cpen_projection, old_space, nmnn_sides);
  end
    
  %----------------------------------------------------------------------
  % update
  time = time + dt;
  u_n = u_n1;
  udot_n = udot_n1;
    
  % check max time
  if (time + dt > problem_data.Time_max)
    dt = problem_data.Time_max - time;
  end

  % Store results
  if (time >= time_save(1))
    save_results_step(u_n1, udot_n1, time, hspace, hmsh, geometry, save_info, save_id)
    results.time(save_id) = time;
    save_id = save_id + 1;
    time_save = time_save(time_save > time);
  end

end % end loop over time steps

disp('----------------------------------------------------------------')
disp(strcat('END ANALYSIS t=',num2str(time)))
disp('----------------------------------------------------------------')

% Crop results
results.time = results.time(1:save_id-1);

end

%--------------------------------------------------------------------------
% FUNCTIONS
%--------------------------------------------------------------------------
% save and plot results
function save_results_step (field, field_dot,time, hspace, hmsh, geometry, save_info, counter)

  if (~isfield(save_info, 'vtk_pts') || isempty (save_info.vtk_pts))
    vtk_pts = {linspace(0, 1, 50), linspace(0, 1, 50)};
  else
    vtk_pts = save_info.vtk_pts;
  end
  
  % save numerical results in a file
  output_file = strcat (save_info.folder_name, '/', save_info.file_name, num2str(counter));
  
  fprintf ('The result is saved in the file %s \n \n', output_file);
  sp_to_vtk (field, hspace, geometry, vtk_pts, output_file ,{'solution', 'gradient'}, {'value', 'gradient'})
  
  save(output_file, 'field', 'field_dot', 'time', 'hspace','hmsh')
  
  fig = figure('visible','off');
  sp_plot_solution (field, hspace, geometry,vtk_pts)
  alpha(0.65)
  shading interp
  hold on
  colorbar
  
  hmsh_plot_cells (hmsh)
  
  title(strcat ('Time: ', num2str(time)))
  view(0,90);
  grid off
  
  saveas(fig , strcat(save_info.folder_name,'/mesh_time_',num2str(time),'.png') )
  close(fig)
end
