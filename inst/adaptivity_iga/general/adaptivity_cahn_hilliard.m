%%%%%%%
%%%%%%%
%%%%%%%
%%%%%%%

function [geometry, hmsh, hspace, results] = ...
              adaptivity_cahn_hilliard (problem_data, method_data, adaptivity_data, initial_conditions, save_info)

%%-------------------------------------------------------------------------
% Initialization of some auxiliary variables
nel = zeros (1, adaptivity_data.num_max_iter); ndof = nel;

% Initialization of the most coarse level of the hierarchical mesh and space
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
% initial conditions
mass_mat = op_u_v_hier (hspace,hspace,hmsh);
% penalty term (matrix)
[Pen, ~] = penalty_matrix (hspace, hmsh, method_data.pen_projection, nmnn_sides);
mass_proj = mass_mat + Pen;

if (isfield(initial_conditions,'fun_u'))
  if (isnumeric(initial_conditions.fun_u) == 1)
    u_0 = fun_u;
  else
    u_0 = compute_initial_conditions(mass_proj, initial_conditions.fun_u, hspace, hmsh);
  end
else
  u_0 = zeros(hspace.ndof,1);
end

if (isfield(initial_conditions,'fun_udot'))
  if (isnumeric(initial_conditions.fun_udot) == 1)
    udot_0 = fun_udot;
  else
    udot_0 = compute_initial_conditions(mass_proj, initial_conditions.fun_udot, hspace, hmsh);
  end
else
  udot_0 = zeros(hspace.ndof,1);
end
clear mass_proj

norm_flux_initial = check_flux_phase_field(hspace, hmsh, u_0, nmnn_sides);
disp(strcat('initial flux =',num2str(norm_flux_initial)))

u_n = u_0;
udot_n = udot_0;

%%-------------------------------------------------------------------------
% generalized alpha parameters
rho_inf_gen_alpha = method_data.rho_inf_gen_alpha;
a_m = .5*((3-rho_inf_gen_alpha)/(1+rho_inf_gen_alpha));
a_f =  1/(1+rho_inf_gen_alpha);
gamma =  .5 + a_m - a_f;

%%-------------------------------------------------------------------------
% loop over time steps

adaptivity_data_flag = true; % if false, the mesh refinement/coarsening is skipped

old_space = struct ('modified', true, 'space', [], 'mesh', [], 'mass_mat', [], ...
  'term3', [], 'term4', [], 'Pen', [], 'pen_rhs', []);

%%-------------------------------------------------------------------------
% Initialize structure to store the results
time = problem_data.initial_time;
save_results_step(u_n, udot_n, time, hspace, hmsh, geometry, save_info, adaptivity_data.estimator_type, 1)
save_id = 2;
results.time = zeros(length(save_info.time_save)+1,1);
flag_stop_save = false;
results.time(1) = time;

dt = method_data.dt;
while time < problem_data.Time_max

  disp('----------------------------------------------------------------')
  disp(strcat('time step t=',num2str(time)))
  disp(strcat('Number of elements = ', num2str(hmsh.nel)))
    
  %----------------------------------------------------------------------
  % adaptivity in space   
  [u_n1, udot_n1, hspace, hmsh, est, old_space] = solve_step_adaptive(u_n, udot_n, hspace, hmsh,  ...
                                              dt, a_m, a_f, gamma, method_data.pen_nitsche, problem_data, ...
                                              adaptivity_data, adaptivity_data_flag, old_space, nmnn_sides);

  %----------------------------------------------------------------------
  % coarsening
  if (time >= adaptivity_data.time_delay)
    if (adaptivity_data_flag)
      [u_n1, udot_n1, hspace, hmsh, old_space] = coarsening_algorithm...
        (est, hmsh, hspace, adaptivity_data, u_n1, udot_n1, method_data.pen_projection, old_space, nmnn_sides);
    end
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
  % update time
  time = time + dt;

  % update
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
%%%%%%%%%%%%%%%%%%%%%%%%%% TODO TODO TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% multiple refinements marking all the elements 
%--------------------------------------------------------------------------
function [hmsh, hspace] = uniform_refinement(hmsh, hspace, n_refinements, adaptivity_data)

  if (n_refinements >= 1)
    for ii = 1:n_refinements
      if (hmsh.nlevels >= adaptivity_data.max_level) % check max depth
        disp('Uniform refinement limited by max_level')
        break
      end

      marked = cell (hmsh.nlevels,1);
      marked{hmsh.nlevels} = hmsh.active{hmsh.nlevels};
      [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
    end
  end
end

%--------------------------------------------------------------------------
% initial conditions
%--------------------------------------------------------------------------
function u_n = compute_initial_conditions(mass_mat, fun_u, hspace, hmsh)
  rhs = op_f_v_hier(hspace,hmsh, fun_u);
  u_n = mass_mat \ rhs;
end

%--------------------------------------------------------------------------
% adaptivity in space
%--------------------------------------------------------------------------
function [u_n1, udot_n1, hspace, hmsh, est, old_space] = solve_step_adaptive...
    (u_n, udot_n, hspace, hmsh, dt, a_m, a_f, gamma, pen, problem_data, ...
    adaptivity_data, adaptivity_data_flag, old_space, nmnn_sides)

  lambda = problem_data.lambda;
  mu = problem_data.mu;
  dmu = problem_data.dmu;
  iter = 0;
  while (1)
    iter = iter + 1;
    disp(strcat('%%%%%%%%%%%%%%%%% Adaptivity iteration ',num2str(iter),' %%%%%%%%%%%%%%%%%'));
    
    %------------------------------------------------------------------
    % solve
    [u_n1, udot_n1, old_space] = generalized_alpha_step(u_n, udot_n, dt, a_m, a_f, gamma, lambda, mu, dmu, ...
                                                        pen, hspace, hmsh, old_space, nmnn_sides);

    %------------------------------------------------------------------
    % skip adaptivity
    if (~adaptivity_data_flag)
      est = zeros(hmsh.nel, 1);
      disp('No adaptivity')
      break
    end

    %------------------------------------------------------------------
    % estimate
%%%%%%%%%%%%%%%%%%%%%%%%%% TODO: write the help of the function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    est = adaptivity_estimate_cahn_hilliard (u_n1, hmsh, hspace, adaptivity_data);

    %------------------------------------------------------------------
    % stopping criteria
    if (iter == adaptivity_data.num_max_iter)
      disp('Warning: reached the maximum number of iterations')
      break
    elseif (hmsh.nlevels > adaptivity_data.max_level)
      disp(strcat('number of levels =',num2str(hmsh.nlevels))) 
      disp('Warning: reached the maximum number of levels')
      break                
    end
        
    if (hmsh.nlevels == adaptivity_data.max_level && hmsh.nel == hmsh.nel_per_level(end))
      disp('Mesh completely refined at the maximum level')
      break
    end

    %------------------------------------------------------------------
    % mark
%%%%%%%%%%%%%%%%%%%%%%%%%% TODO TODO TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    [marked, num_marked_ref] = adaptivity_mark_cahn_hilliard (est, hmsh, hspace, adaptivity_data);
    % limit the maximum refinement depth
    if (hmsh.nlevels == adaptivity_data.max_level)
      num_deleted = numel(marked{hmsh.nlevels});
      marked{hmsh.nlevels} = [];
      num_marked_ref = num_marked_ref - num_deleted ;
    end

    % stopping criterion
    if (num_marked_ref == 0)
      disp('No element is refined')
      break  
    end

    % refine
    [hmsh, hspace, Cref] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
    old_space = struct ('modified', true, 'space', [], 'mesh', [], 'mass_mat', [], ...
                        'term3', [], 'term4', [], 'Pen', [], 'pen_rhs', []);

    % recompute control variables
    [u_n, udot_n] = compute_control_variables_new_mesh(u_n, udot_n, Cref);

  end % end loop adaptivity
end

%--------------------------------------------------------------------------
% coarsening
%--------------------------------------------------------------------------
function [u_n1, udot_n1, hspace, hmsh, old_space] = ...
  coarsening_algorithm(est, hmsh, hspace, adaptivity_data, u_n1, udot_n1, pen_proje, old_space, nmnn_sides)

  %------------------------------------------------------------------  
  % mark
%%%%%%%%%%%%%%%%%%%%%%%%%% TODO TODO TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  [marked, num_marked_coa] = adaptivity_mark_coarsening_cahn_hilliard (est, hmsh, hspace, adaptivity_data);
    
  if (num_marked_coa == 0)
    old_space.modified = false;
  else
  %------------------------------------------------------------------
  % coarsening
    hmsh_fine = hmsh;
    hspace_fine = hspace;
    [hmsh, hspace] = adaptivity_coarsen (hmsh, hspace, marked, adaptivity_data);    
        
    if (hspace.ndof == hspace_fine.ndof)
      old_space.modified = false;
    else
      % recompute control variables:  (mass+penalty) \ (G * u)
      [u_n1, udot_n1] = compute_control_variables_coarse_mesh(hmsh, hspace, hmsh_fine, hspace_fine, u_n1, udot_n1, pen_proje, nmnn_sides);

      old_space = struct ('modified', true, 'space', [], 'mesh', [], 'mass_mat', [], ...
                          'term3', [], 'term4', [], 'Pen', [], 'pen_rhs', []);
    end
    clear hmsh_fine hspace_fine

  end
end

%--------------------------------------------------------------------------
% compute control variables on the new mesh (given the transformation matrix)
%--------------------------------------------------------------------------
function [u_n_ref, udot_n_ref] = compute_control_variables_new_mesh(u_n, udot_n, Cref)
  u_n_ref = Cref * u_n;       
  udot_n_ref = Cref * udot_n;
end

%--------------------------------------------------------------------------
% compute control variables on the coarser mesh by means of L2-projection
%--------------------------------------------------------------------------
function [u_n_coa, udot_n_coa] = ...
   compute_control_variables_coarse_mesh(hmsh, hspace, hmsh_fine, hspace_fine, u_n, udot_n, pen_proje, nmnn_sides)

  mass_coarse = op_u_v_hier(hspace,hspace,hmsh);

  % penalty term (matrix and vector)
  [Pen, ~] = penalty_matrix (hspace, hmsh, pen_proje, nmnn_sides);
  mass_coarse = mass_coarse + Pen;

  hspace_in_hmsh_fine = hspace_in_finer_mesh(hspace, hmsh, hmsh_fine);
  rhs_u = op_Gu_hier (hspace_in_hmsh_fine, hmsh_fine, hspace_fine, u_n);
  u_n_coa = mass_coarse\rhs_u;

  rhs_udot = op_Gu_hier (hspace_in_hmsh_fine, hmsh_fine, hspace_fine, udot_n);
  udot_n_coa = mass_coarse\rhs_udot;

end

%--------------------------------------------------------------------------
% Compute a rhs-vector of the form b_i = (B_i, u) = (B_i, sum_j u_j D_j), 
% with B_i a basis of one (coarse) space, and u written as a combination of 
% basis functions of a second (fine) space, using the finer mesh
%--------------------------------------------------------------------------
function rhs = op_Gu_hier (hspace, hmsh_fine, hspace_fine, uhat_fine)

  rhs = zeros (hspace.ndof, 1);

  ndofs = 0;
  for ilev = 1:hmsh_fine.nlevels
    ndofs = ndofs + hspace.ndof_per_level(ilev);
    if (hmsh_fine.nel_per_level(ilev) > 0)
      x = cell (hmsh_fine.rdim, 1);
      for idim = 1:hmsh_fine.rdim
        x{idim} = reshape (hmsh_fine.msh_lev{ilev}.geo_map(idim,:,:), hmsh_fine.mesh_of_level(ilev).nqn, hmsh_fine.nel_per_level(ilev));
      end

      sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh_fine.msh_lev{ilev}, 'value', true);      
      sp_lev = change_connectivity_localized_Csub (sp_lev, hspace, ilev);

      % coefficients
      ndof_until_lev = sum (hspace_fine.ndof_per_level(1:ilev));
      uhat_lev = hspace_fine.Csub{ilev} * uhat_fine(1:ndof_until_lev);
      sp_lev_fine = sp_evaluate_element_list (hspace_fine.space_of_level(ilev), hmsh_fine.msh_lev{ilev}, 'value', true); 
      sp_lev_fine = change_connectivity_localized_Csub (sp_lev_fine, hspace_fine, ilev);
      utemp = sp_eval_msh (uhat_lev, sp_lev_fine, hmsh_fine.msh_lev{ilev}, {'value'});
      u_fine = utemp{1};

      b_lev = op_f_v (sp_lev, hmsh_fine.msh_lev{ilev}, u_fine);

      dofs = 1:ndofs;
      rhs(dofs) = rhs(dofs) + hspace.Csub{ilev}.' * b_lev;
    end
  end

end

%--------------------------------------------------------------------------
% generalized alpha
%--------------------------------------------------------------------------
function [u_n1, udot_n1, old_space] = generalized_alpha_step(u_n, udot_n, dt, a_m, a_f, gamma, ...
                                       lambda, mu, dmu, pen, hspace, hmsh, old_space, nmnn_sides)

% convergence criteria
  n_max_iter = 20;
  tol_rel_res = 1e-10;
  tol_abs_res = 1e-10;

% predictor step
  u_n1 = u_n;
  udot_n1 = (gamma-1)/gamma * udot_n; 

  for iter = 0:n_max_iter % newton loop

    % field at alpha level
    udot_a = udot_n + a_m *(udot_n1-udot_n);
    u_a = u_n + a_f *(u_n1-u_n);

    % compute the residual (internal)
    [Res_gl, stiff_mat, mass_mat, old_space] = ...
      Res_K_cahn_hilliard(hspace, hmsh, lambda, pen, u_a, udot_a, mu, dmu, old_space, nmnn_sides);

    % convergence check
    if (iter == 0)
      norm_res_0 = norm(Res_gl);
    end
    norm_res =norm(Res_gl);
    
    if (norm_res/norm_res_0 < tol_rel_res) %relative tolerance
      disp(strcat('iteration n°=',num2str(iter)))
      disp(strcat('norm residual=',num2str(norm_res)))
      break
    end
    if (norm_res<tol_abs_res) % absolute tolerance
      disp(strcat('iteration n°=',num2str(iter)))
      disp(strcat('norm residual=',num2str(norm_res)))
      break
    end
    if (iter == n_max_iter) % maximum iterations
        disp(strcat('stop max number of iterations'))
        disp(strcat('norm residual=',num2str(norm_res)))
    end
    
    % Compute the update, and update the solution
    A_gl = a_m * mass_mat + a_f * gamma *dt * stiff_mat ; 
    d_udot = - A_gl\Res_gl;

    udot_n1 = udot_n1 + d_udot;
    u_n1 = u_n1 + gamma * dt* d_udot;
  end

end

%--------------------------------------------------------------------------
% Cahn-Hilliard equation residual and tangent matrix
%--------------------------------------------------------------------------
function [Res_gl, stiff_mat, mass_mat, old_space] = Res_K_cahn_hilliard(hspace, hmsh, lambda, pen, u_a, udot_a, ...
                                                     mu, dmu, old_space, nmnn_sides)

  if (old_space.modified == true)
    
    % mass matrix
    mass_mat = op_u_v_hier (hspace,hspace,hmsh);
    
    % double well (matrices)
    [term2, term2K] = op_gradfu_gradv_hier (hspace, hmsh, u_a, mu, dmu);
     
    % laplacian (matrix)
    term3 = op_laplaceu_laplacev_hier (hspace, hspace, hmsh, lambda);
    
    % term 4 (matrix)
    term4 = int_term_4 (hspace, hmsh, lambda, nmnn_sides);
    
    % penalty term (matrix and vector)
    [Pen, pen_rhs] = penalty_matrix (hspace, hmsh, pen, nmnn_sides);

    % update old_space
    old_space = struct ('modified', false, 'space', hspace, 'mesh', hmsh, ...
      'mass_mat', mass_mat, 'term3', term3, 'term4', term4, 'Pen', Pen, 'pen_rhs', pen_rhs);

  elseif (old_space.modified == false)
    mass_mat =  old_space.mass_mat;
    term3 =  old_space.term3;
    term4 =  old_space.term4;
    Pen = old_space.Pen;
    pen_rhs =  old_space.pen_rhs;

    [term2, term2K] = op_gradfu_gradv_hier(old_space.space, old_space.mesh, u_a, mu, dmu);
  end

  % residual
  Res_gl = mass_mat*udot_a + term2*u_a + term3*u_a - (term4 + term4')*u_a + Pen*u_a - pen_rhs;

  % tangent stiffness matrix (mass is not considered here)
  stiff_mat = term2 + term2K + term3 - (term4 + term4') + Pen;

end

%--------------------------------------------------------------------------
% term 4
%--------------------------------------------------------------------------
function A = int_term_4 (hspace, hmsh, lambda, nmnn_sides)
  A =  spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);    
  for iside = 1:length(nmnn_sides)   
    A = A + op_gradv_n_laplaceu_hier(hspace, hmsh, nmnn_sides(iside), lambda );        
  end
end

%--------------------------------------------------------------------------
% OP_GRADVN_LAPLACEU_HIER: assemble the matrix A = [a(i,j)], a(i,j) = (epsilon (grad v n)_j, Delta u_i), 
%  with n the normal vector to the boundary.
%--------------------------------------------------------------------------
function varargout = op_gradv_n_laplaceu_hier (hspace, hmsh, nmnn_side, coeff)

  if (nargin == 3)
    coeff = @(varargin) ones (size(varargin{1}));
  end
    
  % side with info from the interior for space reconstruction
  hmsh_side_int = hmsh_boundary_side_from_interior (hmsh, nmnn_side);
  shifting_indices = cumsum ([0 hmsh.boundary(nmnn_side).nel_per_level]);
    
  K = spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);

  ndofs_u = 0;
  ndofs_v = 0;
  for ilev = 1:hmsh.boundary(nmnn_side).nlevels

    ndofs_u = ndofs_u + hspace.ndof_per_level(ilev);
    ndofs_v = ndofs_v + hspace.ndof_per_level(ilev);

    if (hmsh.boundary(nmnn_side).nel_per_level(ilev) > 0)
      % mesh of the selected side
      elements = shifting_indices(ilev)+1:shifting_indices(ilev+1); 
      hmsh_side = hmsh_eval_boundary_side (hmsh, nmnn_side, elements);

      x = cell (hmsh.rdim,1);
      for idim = 1:hmsh.rdim
        x{idim} = reshape (hmsh_side.geo_map(idim,:,:), hmsh_side.nqn, hmsh_side.nel);
      end
        
      % reconstruct the part of the space defined on the selected boundary
      msh_side_from_interior = hmsh_side_int.mesh_of_level(ilev);
      sp_bnd = hspace.space_of_level(ilev).constructor (msh_side_from_interior);
      msh_side_from_interior_struct = msh_evaluate_element_list (msh_side_from_interior, hmsh_side_int.active{ilev});
        
      % evaluate the space
      spu_lev = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'laplacian', true);
      spv_lev = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', false, 'gradient', true);
        
      spu_lev = change_connectivity_localized_Csub (spu_lev, hspace, ilev);
      spv_lev = change_connectivity_localized_Csub (spv_lev, hspace, ilev);     

      % compute the matrix
      K_lev = op_gradv_n_laplaceu (spu_lev, spv_lev, hmsh_side, coeff (x{:}));

      dofs_u = 1:ndofs_u;
      dofs_v = 1:ndofs_v;
      Ktmp =  hspace.Csub{ilev}.' * K_lev * hspace.Csub{ilev};
      K(dofs_v,dofs_u) = K(dofs_v,dofs_u) +Ktmp;    
    end
  end

  if (nargout == 1)
    varargout{1} = K;
  elseif (nargout == 3)
    [rows, cols, vals] = find (K);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end

% function varargout = op_gradv_n_laplaceu (spu, spv, msh, coeff)
% 
%   gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], ...
% 		   msh.nqn, spv.nsh_max, msh.nel);
% 
%   ndim = size (gradv, 2);
% 
%   shpu = reshape (spu.shape_function_laplacians, spu.ncomp, msh.nqn, spu.nsh_max, msh.nel);
% 
%   rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
%   cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
%   values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
% 
%   jacdet_weights = msh.jacdet .* msh.quad_weights .* coeff;
% 
%   ncounter = 0;
%   for iel = 1:msh.nel
%     if (all (msh.jacdet(:,iel)))
%       gradv_iel = gradv(:,:,:,:,iel);
%       normal_iel = reshape (msh.normal(:,:,iel), [1, ndim, msh.nqn]);
% 
%       gradv_n = reshape (sum (bsxfun (@times, gradv_iel, normal_iel), 2), spv.ncomp, msh.nqn, spv.nsh_max, 1);
%       shpu_iel = reshape (shpu(:, :, :, iel), spu.ncomp, msh.nqn, 1, spu.nsh_max);
% 
%       jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);
% 
%       gradv_n_times_jw = bsxfun (@times, jacdet_iel, gradv_n);
%       tmp1 = sum (bsxfun (@times, gradv_n_times_jw, shpu_iel), 1);
%       elementary_values = reshape (sum (tmp1, 2), spv.nsh_max, spu.nsh_max);
% 
%       [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spu.connectivity(:,iel));
%       indices = rows_loc & cols_loc;
%       rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc(indices);
%       cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc(indices);
%       values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = elementary_values(indices);
%       ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);
% 
%     else
%       warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradv_n_u: singular map in element number %d', iel)
%     end
%   end
% 
%   if (nargout == 1 || nargout == 0)
%     varargout{1} = sparse (rows, cols, values, spv.ndof, spu.ndof);
%   elseif (nargout == 3)
%     varargout{1} = rows;
%     varargout{2} = cols;
%     varargout{3} = values;
%   else
%     error ('op_gradv_n_u: wrong number of output arguments')
%   end
% 
% end

% %--------------------------------------------------------------------------
% % integral of the double-well function
% %--------------------------------------------------------------------------
% function [A, B] = op_gradfu_gradv_hier (hspace, hmsh, uhat, f, df)
% 
%   A = spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);
%   B = spalloc (hspace.ndof, hspace.ndof, 6*hspace.ndof);
% 
%   ndofs_u = 0;
%   ndofs_v = 0;
%   for ilev = 1:hmsh.nlevels
%     ndofs_u = ndofs_u + hspace.ndof_per_level(ilev);
%     ndofs_v = ndofs_v + hspace.ndof_per_level(ilev);
%     if (hmsh.nel_per_level(ilev) > 0)
% 
%       % space
%       spu_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true, 'gradient', true, 'laplacian', true);
%       spu_lev = change_connectivity_localized_Csub (spu_lev, hspace, ilev);
% 
%       % coefficients
%       ndof_until_lev = sum (hspace.ndof_per_level(1:ilev));
%       uhat_lev = hspace.Csub{ilev} * uhat(1:ndof_until_lev);
% 
%       utemp = sp_eval_msh (uhat_lev, spu_lev, hmsh.msh_lev{ilev}, {'value', 'gradient'});
%       u = utemp{1};
%       gradu = utemp{2};
% 
%       coeffs_A = f(u);
%       coeffs_B = df(u);
%       coeffs_Bv = gradu;
%       for idim = 1:hmsh.ndim
%         coeffs_Bv(idim,:,:)=coeffs_Bv(idim,:,:) .* reshape(coeffs_B, 1, size(coeffs_B,1), size(coeffs_B,2) );
%       end
% 
%       % matrices
%       A_lev = op_gradu_gradv (spu_lev, spu_lev, hmsh.msh_lev{ilev}, coeffs_A);
%       B_lev = op_vel_dot_gradu_v (spu_lev, spu_lev, hmsh.msh_lev{ilev}, coeffs_Bv)';
% 
%       % assemble
%       dofs_u = 1:ndofs_u;
%       dofs_v = 1:ndofs_v;
% 
%       A(dofs_v,dofs_u) = A(dofs_v,dofs_u) + hspace.Csub{ilev}.' * A_lev * hspace.Csub{ilev};   
%       B(dofs_v,dofs_u) = B(dofs_v,dofs_u) + hspace.Csub{ilev}.' * B_lev * hspace.Csub{ilev};   
%     end
%   end
% 
% end

% %--------------------------------------------------------------------------
% % operator for boundary conditions
% %--------------------------------------------------------------------------
% function varargout = op_gradv_n_u_hier (hspace, hmsh, nmnn_side, coeff)
% 
%   if (nargin == 3)
%     coeff = @(varargin) ones (size(varargin{1}));
%   end
% 
%   % side with info from the interior for space reconstruction
%   hmsh_side_int = hmsh_boundary_side_from_interior (hmsh, nmnn_side);
%   shifting_indices = cumsum ([0 hmsh.boundary(nmnn_side).nel_per_level]);
% 
%   K = spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);
% 
%   ndofs_u = 0;
%   ndofs_v = 0;
%   for ilev = 1:hmsh.boundary(nmnn_side).nlevels
% 
%     ndofs_u = ndofs_u + hspace.ndof_per_level(ilev);
%     ndofs_v = ndofs_v + hspace.ndof_per_level(ilev);
% 
%     if (hmsh.boundary(nmnn_side).nel_per_level(ilev) > 0)
% 
%       % mesh of the selected side
%       elements = shifting_indices(ilev)+1:shifting_indices(ilev+1); 
%       hmsh_side = hmsh_eval_boundary_side (hmsh, nmnn_side, elements);
% 
%       x = cell (hmsh.rdim,1);
%       for idim = 1:hmsh.rdim
%         x{idim} = reshape (hmsh_side.geo_map(idim,:,:), hmsh_side.nqn, hmsh_side.nel);
%       end
% 
%       % reconstruct the part of the space defined on the selected boundary
%       msh_side_from_interior = hmsh_side_int.mesh_of_level(ilev);
%       sp_bnd = hspace.space_of_level(ilev).constructor (msh_side_from_interior);
%       msh_side_from_interior_struct = msh_evaluate_element_list (msh_side_from_interior, hmsh_side_int.active{ilev});
% 
%       % evaluate the space
%       spu_lev = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', true);
%       spv_lev = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', false, 'gradient', true);
% 
%       spu_lev = change_connectivity_localized_Csub (spu_lev, hspace, ilev);
%       spv_lev = change_connectivity_localized_Csub (spv_lev, hspace, ilev);     
% 
%       % compute the matrix
%       K_lev = op_gradv_n_u (spu_lev, spv_lev, hmsh_side, coeff (x{:}));
% 
%       dofs_u = 1:ndofs_u;  
%       dofs_v = 1:ndofs_v;
%       Ktmp =  hspace.Csub{ilev}.' * K_lev * hspace.Csub{ilev};
%       K(dofs_v,dofs_u) = K(dofs_v,dofs_u) +Ktmp;    
%     end
%   end
% 
%   if (nargout == 1)
%     varargout{1} = K;
%   elseif (nargout == 3)
%     [rows, cols, vals] = find (K);
%     varargout{1} = rows;
%     varargout{2} = cols;
%     varargout{3} = vals;
%   end
% 
% end

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

%---------------------------

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
% penalty term
%--------------------------------------------------------------------------
function [P, rhs] = penalty_matrix (hspace, hmsh, pen, nmnn_sides)
    
  P =  spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);
  rhs = zeros(hspace.ndof,1);

  for iside = 1:length(nmnn_sides)
    [mass_pen, rhs_pen] = penalty_grad (hspace, hmsh, nmnn_sides(iside), pen);
    P = P + mass_pen; 
    rhs = rhs + rhs_pen;
  end

end

%--------------------------------------------------------------------------
function [mass_pen, rhs_pen] = penalty_grad (hspace, hmsh, nmnn_side, pen)

  hmsh_side_int = hmsh_boundary_side_from_interior (hmsh, nmnn_side ) ;
  shifting_indices = cumsum ([0 hmsh.boundary(nmnn_side).nel_per_level]);

  mass_pen = spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);

  ndofs_u = 0;
  ndofs_v = 0;
  for ilev = 1:hmsh.boundary(nmnn_side).nlevels

    ndofs_u = ndofs_u + hspace.ndof_per_level(ilev);
    ndofs_v = ndofs_v + hspace.ndof_per_level(ilev);
    
    if (hmsh.boundary(nmnn_side).nel_per_level(ilev) > 0)

      % mesh of the selected side
      elements = shifting_indices(ilev)+1:shifting_indices(ilev+1); 
      hmsh_side = hmsh_eval_boundary_side (hmsh, nmnn_side, elements);

      % reconstruct the part of the space defined on the selected boundary
      msh_side_from_interior = hmsh_side_int.mesh_of_level(ilev);
      sp_bnd = hspace.space_of_level(ilev).constructor (msh_side_from_interior);
      msh_side_from_interior_struct = msh_evaluate_element_list (msh_side_from_interior, hmsh_side_int.active{ilev});

      % evaluate the space
      spu_lev = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'gradient', true);
      spu_lev = change_connectivity_localized_Csub (spu_lev, hspace, ilev);

      % compute the matrix
      coe_side = repmat(hmsh_side.element_size, hmsh_side.nqn,1);
      coe_side = pen .* coe_side;
      K_lev = op_gradu_n_gradv_n(spu_lev, spu_lev, hmsh_side, coe_side);

      dofs_u = 1:ndofs_u;
      dofs_v = 1:ndofs_v;
      Ktmp =  hspace.Csub{ilev}.' * K_lev * hspace.Csub{ilev};
      mass_pen(dofs_v,dofs_u) = mass_pen(dofs_v,dofs_u) + Ktmp;
    end
  end

  rhs_pen  = zeros(hspace.ndof,1);

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


% %--------------------------------------------------------------------------
% % TODO: WHAT DOES THIS FUNCTION DO?
% %--------------------------------------------------------------------------
% 
% function [adaptivity_data_flag, last_mesh] = evaluate_adaptivity_likelihood(hmsh, hspace, u, last_mesh, adaptivity_data)
% 
% u_old = last_mesh.u;
% 
% level_shared = adaptivity_data.max_level-1;
% [hmsh, hspace, u] = refine_upto_level(hmsh, hspace, u, level_shared, adaptivity_data );
% 
% uex     = @(x, y) zeros (size (x));
% graduex = @(x, y) cat (1, ...
%                    reshape ( zeros (size (x)) , [1, size(x)]), ...
%                    reshape ( zeros (size (x)) , [1, size(x)]));
% 
% [~, ~, nume, ~, ~, ~] = sp_h1_error (hspace, hmsh, u-u_old, uex, graduex);
% [~, ~, deno, ~, ~, ~] = sp_h1_error (hspace, hmsh, u, uex, graduex);
% err = nume/deno;
% 
% if (err> last_mesh.tol_skip_ada)
%   adaptivity_data_flag = true;
%   last_mesh.mesh = hmsh;
%   last_mesh.space = hspace;
%   last_mesh.u = u;
% else
%   adaptivity_data_flag = false;
% end
% 
% end
% 
% 
% function [hmsh, hspace, u] = refine_upto_level(hmsh, hspace, u, level, adaptivity_data)
% 
%   for ilev = 1:level
%     marked = cell(hmsh.nlevels,1); 
% 
%     elements = hmsh.active{ilev};  
%     marked{ilev} = elements;
% 
%     [hmsh, hspace, Cref] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
%     u = Cref*u;
%   end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [hmsh, hspace, Ccoar] = my_adaptivity_coarsen (hmsh, hspace, marked, adaptivity_data)
% 
% switch (adaptivity_data.flag)
%   case 'functions'
%     error ('Coarsening by functions is not implemented yet')
%   case 'elements'
%     marked_elements = marked;
% end
% [reactivated_elements, ~] = mark_elements_to_reactivate_from_active (marked_elements, hmsh, hspace, adaptivity_data);
% 
% hmsh_fine = hmsh;
% [hmsh, removed_cells] = hmsh_coarsen (hmsh, reactivated_elements);
% 
% reactivated_fun = functions_to_reactivate_from_cells (hmsh, hspace, reactivated_elements);
% 
% if (nargout == 3)
%   hspace_fine = hspace;
%   hspace = hspace_coarsen (hspace, hmsh, reactivated_fun, removed_cells);
%   warning ('Coarsening matrix computed with an expensive L^2 projection')
%   M = op_u_v_hier (hspace, hspace, hmsh);
%   G = op_u_v_hier (hspace_fine, hspace_in_finer_mesh(hspace, hmsh, hmsh_fine), hmsh_fine);
%   Ccoar = M \ G; Ccoar(abs(Ccoar) < 1e-12) = 0;
% else
%   hspace = hspace_coarsen (hspace, hmsh, reactivated_fun, removed_cells);
% end
% 
% hmsh = hmsh_remove_empty_levels (hmsh);
% hspace = hspace_remove_empty_levels (hspace, hmsh);
% 
% end



% function [deact_marked, num] = mark_elements_to_reactivate_from_active (marked, hmsh, hspace, adaptivity_data)
% 
% if (~isfield (adaptivity_data, 'coarsening_flag'))
%   adaptivity_data.coarsening_flag = 'all';
% end
% 
% if (~isfield(adaptivity_data,'adm_class') || adaptivity_data.adm_class < 2)
%   adm_class = 0;
% else
%   adm_class = adaptivity_data.adm_class;
% end
% 
% if (~isfield(adaptivity_data,'adm_type') || isempty (adaptivity_data.adm_type))
%   adm_type = 'T-admissible';
% else
%   adm_type = adaptivity_data.adm_type;
% end
% 
% 
% deact_marked = cell (hmsh.nlevels, 1);
% 
% for lev = hmsh.nlevels-1:-1:1
%   if (~isempty(marked{lev+1}))
%     [parents, flag] = hmsh_get_parent (hmsh, lev+1, marked{lev+1});
%     if (flag ~= 1)
%       error ('Some nonactive elements were marked.')
%     end
% 
%     [~,~,children_per_cell] = hmsh_get_children (hmsh, lev, parents);
% 
%     ind = all (ismember (children_per_cell, hmsh.active{lev+1}));
%     parents = parents(ind);
%     children_per_cell = children_per_cell(:,ind);
% 
%     if (strcmpi (adaptivity_data.coarsening_flag, 'any'))
%       deact_marked{lev} = parents;
%     elseif (strcmpi (adaptivity_data.coarsening_flag, 'all'))
%       ind2 = all (ismember (children_per_cell, marked{lev+1}));
%       parents = parents(ind2);
%       children_per_cell = children_per_cell(:,ind2);
%       deact_marked{lev} = parents;
%     else
%       error ('Unknown option for coarsening, in adaptivity_data.coarsening_flag')
%     end
%     marked{lev+1} = children_per_cell;
% 
% % Algorithm to recover admissible meshes
%     if (adm_class)
%       lev_s = lev + adm_class;
%       if (lev_s > hmsh.nlevels)
%         continue
%       else
%         active_and_deact = union (hmsh.active{lev_s}, hmsh.deactivated{lev_s});
%         active_and_deact = setdiff (active_and_deact, marked{lev_s});
%         keep_inds = [];
% 
%         if (strcmpi (adm_type, 'T-admissible'))
%           supp_ext = support_extension (hmsh, hspace, children_per_cell(:), lev+1, lev+1);
%           [~, descendants_of_cell] = hmsh_get_descendants (hmsh, supp_ext, lev+1, lev_s);
%           for iel = 1:numel(deact_marked{lev})
%             supp_ext_local = support_extension (hmsh, hspace, children_per_cell(:,iel), lev+1, lev+1);
%             [~,ia,~] = intersect (supp_ext, supp_ext_local);
%             if (isempty (intersect (descendants_of_cell(:,ia), active_and_deact)))
%               keep_inds = [keep_inds, iel];
%             end          
%           end
% 
%         elseif (strcmpi (adm_type, 'H-admissible'))
%           supp_ext = support_extension (hmsh, hspace, deact_marked{lev}, lev, lev);
%           [~, descendants_of_cell] = hmsh_get_descendants (hmsh, supp_ext, lev, lev_s);
%           for iel = 1:numel(deact_marked{lev})
%             supp_ext_local = support_extension (hmsh, hspace, deact_marked{lev}(iel), lev, lev);
%             [~,ia,~] = intersect (supp_ext, supp_ext_local);
%             if (isempty (intersect (descendants_of_cell(:,ia), active_and_deact)))
%               keep_inds = [keep_inds, iel];
%             end          
%           end
%         end
%         deact_marked{lev} = deact_marked{lev}(keep_inds);
%       end
%     end
% 
%   end
% end
% 
% num = sum (cellfun (@numel, deact_marked));
% 
% end
