% MP_SOLVE_CAHN_HILLIARD_C1: solve the Cahn-Hilliard equation, with C^1 multipatch splines, and a generalized alpha discretization in time.
%
% The functions solves the problem of finding u such that
%
%  du/dt - Delta (mu(u) - lambda*Delta u) = 0
%
% with Delta the Laplacian, and mu(u) = alpha u^3 - beta u, and Neumann boundary conditions.
%
% The values of alpha and beta (or mu itself) can be changed in op_gradmu_gradv_tp.
%
% For details on the problem and the formulation, see
%  H. Gomez, V.M. Calo, Y. Bazilevs, T.J.R. Hughes, CMAME 197 (2008), 4333-4352.
%  H. Gomez, A. Reali, G. Sangalli, J. Comput. Physics 262 (2014), 153-171.
%
% USAGE:
%
%   [geometry, msh, space, results] = mp_solve_cahn_hilliard_C1 (problem_data, method_data, save_info)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - periodic_directions: parametric directions along which to apply periodic conditions (may be empty)
%    - lambda:       parameter representing the length scale of the problem, and the width of the interface
%    - Time_max:     final time
%    - fun_u:        initial condition. Equal to zero by default.
%    - fun_udot:     initial condition for time derivative. Equal to zero by default.
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions (at most degree minus two).
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%    - dt:         time step size for generalized-alpha method
%    - rho_inf_gen_alpha: parameter in [0,1], which governs numerical damping of the generalized alpha method
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space:    space object that defines the discrete space (see sp_scalar)
%  results:  a struct with the saved results, containing the following fields:
%    - time: (array of length Ntime) time at which the solution was saved
%    - u:    (size ndof x Ntime) degrees of freedom for the solution
%    - udot: (size ndof x Ntime) degrees of freedom for the time derivative
%
% Only periodic and Neumann boundary conditions are implemented. Neumann
%  conditions are considered by default.
%
% Copyright (C) 2023 Michele Torre, Rafael Vazquez
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

function [geometry, hmsh, hspace, results] = adaptivity_cahn_hilliard_mp_C1 (problem_data, method_data, adaptivity_data, initial_conditions, save_info)

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

if (initial_conditions.restart_flag == false)
  % Refine the mesh up to a predefined level
  n_refinements = adaptivity_data.max_level - 1; % number of uniform refinements
  [hmsh, hspace] = uniform_refinement(hmsh, hspace, n_refinements, adaptivity_data);
end

%%-------------------------------------------------------------------------
% Initial conditions, with a penalty term
if (initial_conditions.restart_flag == true)
  disp('restart analysis')
  hspace = initial_conditions.space_reload;
  hmsh = initial_conditions.mesh_reload;
  u_n = initial_conditions.fun_u;
  udot_n = initial_conditions.fun_udot;
  time = initial_conditions.time;
else
  mass_mat = op_u_v_hier (hspace,hspace,hmsh);
  [Pen, ~] = op_penalty_dudn (hspace, hmsh, nmnn_sides, method_data.Cpen_projection);
  mass_proj = mass_mat + Pen;

  if (isfield(initial_conditions,'fun_u') && ~isempty(initial_conditions.fun_u))
    rhs = op_f_v_hier (hspace, hmsh, initial_conditions.fun_u);
    u_n = mass_proj \ rhs;
  else
    u_n = zeros(hspace.ndof, 1);
  end
    
  if (isfield(initial_conditions,'fun_udot') && ~isempty(initial_conditions.fun_udot))
    rhs = op_f_v_hier (hspace, hmsh, initial_conditions.fun_udot);
    udot_n = mass_proj \ rhs;
  else
    udot_n = zeros(hspace.ndof, 1);
  end
  clear mass_proj
end

%%-------------------------------------------------------------------------
% Generalized alpha parameters
rho_inf_gen_alpha = method_data.rho_inf_gen_alpha;
a_m = .5*((3-rho_inf_gen_alpha)/(1+rho_inf_gen_alpha));
a_f =  1/(1+rho_inf_gen_alpha);
gamma =  .5 + a_m - a_f;

%%-------------------------------------------------------------------------
% save matrices previous mesh
old_space = struct ('modified', true, 'space', [], 'mesh', [], 'mass_mat', [], ...
  'lapl_mat', [], 'bnd_mat', [], 'Pen', [], 'pen_rhs', []);

%%-------------------------------------------------------------------------
% Initialize structure to store the results
time = problem_data.initial_time;
save_results_step(u_n, udot_n, time, hspace, hmsh, geometry, save_info, adaptivity_data.estimator_type, 1)
save_id = 2;
results.time = zeros(length(save_info.time_save)+1,1);
flag_stop_save = false;
results.time(1) = time;

%%-------------------------------------------------------------------------
% Loop over time steps
dt = method_data.dt;
while time < problem_data.Time_max

  disp('----------------------------------------------------------------')
  disp(strcat('time step t=',num2str(time)))
  disp(strcat('Number of elements = ', num2str(hmsh.nel)))

  %----------------------------------------------------------------------
  % adaptivity in space
  [u_n1, udot_n1, hspace, hmsh, est, old_space] = solve_step_adaptive_CH(u_n, udot_n, hspace, hmsh,  ...
                                              dt, a_m, a_f, gamma, method_data.Cpen_nitsche, problem_data, ...
                                              adaptivity_data, old_space, nmnn_sides);
    
  %----------------------------------------------------------------------
  % coarsening
  if (time >= adaptivity_data.time_delay)
    [u_n1, udot_n1, hspace, hmsh, old_space] = coarsening_algorithm ...
      (est, hmsh, hspace, adaptivity_data, u_n1, udot_n1, method_data.Cpen_projection, old_space, nmnn_sides);
  end
    
  % Store results
  if (flag_stop_save == false)
    if (time +dt >= save_info.time_save(save_id))  
      save_results_step(u_n1, udot_n1,time+dt, hspace, hmsh, geometry, save_info, adaptivity_data.estimator_type, save_id)
      if (save_id > length(save_info.time_save))
        flag_stop_save = true;
      end
      save_id = save_id + 1;
    end
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
%--------------------------------------------------------------------------
% save and plot results
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_results_step(field, field_dot,time, hspace, hmsh, geometry, save_info, estimator_type, counter)

    
    vtk_pts = {linspace(0, 1, 32*4), linspace(0, 1, 32*4)};
        
    % save numerical results in a file     
        
    output_file = strcat( save_info.folder_name,'/Square_cahn_hilliard_adaptive_',estimator_type,'_', num2str(counter)  );
    fprintf ('The result is saved in the file %s \n \n', output_file);
    sp_to_vtk (field, hspace, geometry, vtk_pts, output_file ,{'solution', 'gradient'}, {'value', 'gradient'})    
    
    save(output_file,'field', 'field_dot', 'time','hspace','hmsh')

%     
%     fig = figure('visible','off');
%     sp_plot_solution (field, hspace, geometry,vtk_pts)
%     alpha(0.5)
%     shading interp
%     hold on
%     colorbar
%     
%     hmsh_plot_cells (hmsh)
%     
%     title(time)
%     view(0,90);
%     grid off
%     saveas(fig , strcat(save_info.folder_name,'/mesh_time_',num2str(time),'.png') )
%     close(fig)

end
