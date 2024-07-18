%%%%%%%
%%%%%%%
%%%%%%%
%%%%%%%

function [geometry, hmsh, hspace, results] = ...
              adaptivity_cahn_hilliard (problem_data, method_data, adaptivity_data, initial_conditions, save_info)

%%-------------------------------------------------------------------------
% Initialization of the coarsest level of the hierarchical mesh and space
[hmsh, hspace, geometry] = adaptivity_initialize_laplace (problem_data, method_data);

if (exist ('nmnn_sides','var') && ~isempty (nmnn_sides))
  disp('User defined Neumann sides deleted. Neumann conditions used everywhere.')
  clear nmnn_sides
end
nmnn_sides = 1:2*hmsh.ndim;

% Refine the mesh up to a predefined level
n_refinements = adaptivity_data.max_level - 1; % number of uniform refinements
[hmsh, hspace] = uniform_refinement(hmsh, hspace, n_refinements, adaptivity_data);

%%-------------------------------------------------------------------------
% initial conditions, with a penalty term
mass_mat = op_u_v_hier (hspace,hspace,hmsh);
[Pen, ~] = penalty_matrix (hspace, hmsh, nmnn_sides, method_data.Cpen_projection);
[Pen2, ~] = op_penalty_dudn (hspace, hmsh, nmnn_sides, method_data.Cpen_projection);
[Pen2, rhs2] = op_penalty_dudn (hspace, hmsh, nmnn_sides, method_data.Cpen_projection, @(x,y) cos(x)+3*sin(y));
mass_proj = mass_mat + Pen;

if (isfield(initial_conditions,'fun_u'))
  if (isnumeric(initial_conditions.fun_u))
    u_n = fun_u;
  else
    rhs = op_f_v_hier(hspace,hmsh, initial_conditions.fun_u);
    u_n = mass_proj \ rhs;
  end
else
  u_n = zeros(hspace.ndof,1);
end

if (isfield(initial_conditions,'fun_udot'))
  if (isnumeric(initial_conditions.fun_udot))
    udot_n = fun_udot;
  else
    rhs = op_f_v_hier(hspace,hmsh, initial_conditions.fun_udot);
    udot_n = mass_proj \ rhs;
  end
else
  udot_n = zeros(hspace.ndof,1);
end
clear mass_proj

% norm_flux_initial = check_flux_phase_field(hspace, hmsh, u_n, nmnn_sides);
% disp(strcat('initial flux =',num2str(norm_flux_initial)))

%%-------------------------------------------------------------------------
% generalized alpha parameters
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
  [u_n1, udot_n1, hspace, hmsh, est, old_space] = solve_step_adaptive(u_n, udot_n, hspace, hmsh,  ...
                                              dt, a_m, a_f, gamma, method_data.Cpen_nitsche, problem_data, ...
                                              adaptivity_data, old_space, nmnn_sides);

  %----------------------------------------------------------------------
  % coarsening
  if (time >= adaptivity_data.time_delay)
    [u_n1, udot_n1, hspace, hmsh, old_space] = coarsening_algorithm...
      (est, hmsh, hspace, adaptivity_data, u_n1, udot_n1, method_data.Cpen_projection, old_space, nmnn_sides);
  end

  % Store results
  if (~flag_stop_save)
    if (time + dt >= save_info.time_save(save_id))
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
    dt = problem_data.Time_max-time;
  end
    
end % end loop over time steps

disp('----------------------------------------------------------------')
disp(strcat('END ANALYSIS t=',num2str(time)))
disp('----------------------------------------------------------------')

% crop results
results.time = results.time(1:save_id-1);

end

%--------------------------------------------------------------------------
% FUNCTIONS
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save and plot results 
function save_results_step (field, field_dot,time, hspace, hmsh, geometry, save_info, estimator_type, counter)


vtk_pts = {linspace(0, 1, 32*4), linspace(0, 1, 32*4)};
    
% save numerical results in a file
      
output_file = strcat( save_info.folder_name,'/Square_cahn_hilliard_adaptive_',estimator_type,'_', num2str(counter)  );
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (field, hspace, geometry, vtk_pts, output_file ,{'solution', 'gradient'}, {'value', 'gradient'})    

save(output_file,'field', 'field_dot', 'time','hspace','hmsh')

fig = figure('visible','off');
% [eu, F] = sp_eval (u, hspace, geometry, npts,{'gradient'}); % sp_eval (u, hspace, geometry, npts,{'value', 'gradient'});
% fun_to_plot = reshape(sqrt(eu(1,:,:).^2 + eu(2,:,:).^2), size(F,2), size(F,3));
[eu, F] = sp_eval (field, hspace, geometry, vtk_pts,{'value'}); % sp_eval (u, hspace, geometry, npts,{'value', 'gradient'});
fun_to_plot = eu;
surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), fun_to_plot, 'FaceAlpha',0.5);
shading interp
hold on
colorbar

hmsh_plot_cells (hmsh)

title(time)
view(0,90);
grid off

saveas(fig , strcat(save_info.folder_name,'/mesh_time_',num2str(time),'.png') )
close(fig)

end

%--------------------------------------------------------------------------
% check flux through the boundaries
%--------------------------------------------------------------------------
function norm_flux = check_flux_phase_field(hspace, hmsh, uhat, sides)
  norm_flux = 0;
  for iside=1:length(sides) 
    norm_flux_side = flux_side(hspace, hmsh, uhat, sides(iside));
    norm_flux = norm_flux + norm_flux_side;
  end
end

%--------------------------------------------------------------------------
function norm_flux_side = flux_side(hspace, hmsh, uhat, side)

  norm_flux_side= 0;

  % side with info from the interior for space reconstruction
  hmsh_side_int = hmsh_boundary_side_from_interior (hmsh, side);
  shifting_indices = cumsum ([0 hmsh.boundary(side).nel_per_level]);

  for ilev = 1:hmsh.boundary(side).nlevels
    if (hmsh.boundary(side).nel_per_level(ilev) > 0)

      % mesh of the selected side
      elements = shifting_indices(ilev)+1:shifting_indices(ilev+1); 
      hmsh_side = hmsh_eval_boundary_side (hmsh, side, elements);

      % reconstruct the part of the space defined on the selected boundary
      msh_side_from_interior = hmsh_side_int.mesh_of_level(ilev);
      sp_bnd = hspace.space_of_level(ilev).constructor (msh_side_from_interior);
      msh_side_from_interior_struct = msh_evaluate_element_list (msh_side_from_interior, hmsh_side_int.active{ilev});
            
      % evaluate the space
      sp_lev = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', false, 'gradient', true);
      sp_lev = change_connectivity_localized_Csub (sp_lev, hspace, ilev);

      ndof_until_lev = sum (hspace.ndof_per_level(1:ilev));
      uhat_lev = hspace.Csub{ilev} * uhat(1:ndof_until_lev);

      gradu = sp_eval_msh (uhat_lev, sp_lev, msh_side_from_interior_struct, 'gradient');
      valu = reshape (sum (gradu .* hmsh_side.normal, 1), hmsh_side.nqn, hmsh_side.nel);

      w = hmsh_side.quad_weights .* hmsh_side.jacdet;
      errl2_elem = sum(valu.*w,1);
      errl2  = sum (errl2_elem);

      norm_flux_side = norm_flux_side + errl2;
    end
  end

end

%--------------------------------------------------------------------------
% plot results
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%% TODO TODO TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_results(geometry, hmsh, hspace, u, time, filenumber)

npts = [51 51];

fig = figure('visible','off');
% [eu, F] = sp_eval (u, hspace, geometry, npts,{'gradient'}); % sp_eval (u, hspace, geometry, npts,{'value', 'gradient'});
% fun_to_plot = reshape(sqrt(eu(1,:,:).^2 + eu(2,:,:).^2), size(F,2), size(F,3));
[eu, F] = sp_eval (u, hspace, geometry, npts,{'value'}); % sp_eval (u, hspace, geometry, npts,{'value', 'gradient'});
fun_to_plot = eu;
surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), fun_to_plot, 'FaceAlpha',0.5);
shading interp
hold on
colorbar

hmsh_plot_cells (hmsh)

title(time)
view(0,90);
grid off

saveas(fig , strcat('results_adaptive/mesh_time_',num2str(time),'.png') )
close(fig)

output_file =  strcat('results_adaptive/cahn_hilliard_adaptive_order_',num2str(filenumber) ) ;
sp_to_vtk (u, hspace, geometry, npts, output_file,{'u', 'grad_u'}, {'value', 'gradient'})

end
