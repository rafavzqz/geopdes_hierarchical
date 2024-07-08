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
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

data_names = fieldnames (initial_conditions);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= initial_conditions.(data_names{iopt});']);
end


%%-------------------------------------------------------------------------
% Parameters for the double-well function (to be given in problem_data)
alpha = 1; 
beta = 1;
mu = @(x) 3 * alpha * x.^2 - beta;
dmu = @(x) 6 * alpha * x;

%%-------------------------------------------------------------------------
% Generalized alpha parameters
a_m = .5*((3-rho_inf_gen_alpha)/(1+rho_inf_gen_alpha));
a_f =  1/(1+rho_inf_gen_alpha);
gamma =  .5 + a_m - a_f;


%%-------------------------------------------------------------------------
% Construct geometry structure, and information for interfaces and boundaries

% Initialization of the most coarse level of the hierarchical mesh and space
[hmsh, hspace, geometry] = adaptivity_initialize_laplace_mp_C1 (problem_data, method_data);

if  initial_conditions.restart_flag == 1
    % skip refinement

else   
    
    % Refine the mesh up to a predefined level
    n_refininements = adaptivity_data.max_level-1; % number of uniform refinements
    [hmsh, hspace] = uniform_refinement(hmsh, hspace, n_refininements, adaptivity_data);

end

%%-------------------------------------------------------------------------
% Initial conditions
time = problem_data.time;

if initial_conditions.restart_flag == 1
    disp('restart analysis')
    hspace = initial_conditions.space_reload;
    hmsh = initial_conditions.mesh_reload;
    u_n = initial_conditions.fun_u;
    udot_n = initial_conditions.fun_udot;
    time = initial_conditions.time;

else

    % mass matrix
    mass_mat = op_u_v_hier(hspace,hspace,hmsh);
    % penalty matrix
    bou   = 1:numel(hmsh.mesh_of_level(1).boundaries);
    [Pen, ~] = penalty_matrix (hspace, hmsh, bou, Cpen_projection);  

    
    if (exist('fun_u', 'var') && ~isempty(fun_u))
      rhs = op_f_v_hier (hspace, hmsh, fun_u);
      u_n = (mass_mat +  Pen)\rhs;
    else
      u_n = zeros(hspace.ndof, 1);
    end
    
    if (exist('fun_udot', 'var') && ~isempty(fun_udot))
      rhs = op_f_v_hier (space, msh, fun_udot);
      udot_n = (mass_mat +  Pen)\rhs;
    else
      udot_n = zeros(hspace.ndof, 1);
    end

end

%%-------------------------------------------------------------------------
% save matrices previous mesh

adaptivity_data_flag = true; % if false, the mesh refinement/coarsening is skipped


old_space.modified = 1;
old_space.space =    [];
old_space.mesh =     [];
old_space.mass_mat = [];
old_space.lapl_mat = [];
old_space.bnd_mat =  [];
old_space.pen   =    [];
old_space.pen_rhs =  [];


counter_plot = 1;


%%-------------------------------------------------------------------------
% Initialize structure to store the results
save_results_step(u_n, udot_n,time, hspace, hmsh, geometry, save_info, adaptivity_data.estimator_type, 1)
save_id = 2;
results.time = zeros(length(save_info.time_save)+1,1);
flag_stop_save = false;
results.time(1) = time;



%%-------------------------------------------------------------------------
% Loop over time steps
while time < Time_max
  disp('----------------------------------------------------------------')
  disp(strcat('time step t=',num2str(time)))

    %----------------------------------------------------------------------
    % adaptivity in space

    [u_n1, udot_n1, hspace, hmsh, est, old_space] = solve_step_adaptive(u_n, udot_n, hspace, hmsh,  ...
                                                dt, a_m, a_f, gamma, lambda, mu, dmu, Cpen_nitsche, ...
                                                problem_data, adaptivity_data, adaptivity_data_flag, old_space);
    
    %----------------------------------------------------------------------
    % coarsening
    if time >= adaptivity_data.time_delay
        if adaptivity_data_flag == true
            [u_n1, udot_n1, hspace, hmsh, old_space] = coarsening_algorithm(est, hmsh, hspace, adaptivity_data,  u_n1, udot_n1, Cpen_projection, old_space);
        end
    end
    
    %----------------------------------------------------------------------
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
    if time + dt >Time_max
        dt = Time_max-time;
    end
    


end % end loop over time steps

disp('----------------------------------------------------------------')
disp(strcat('END ANALYSIS t=',num2str(time)))
disp('----------------------------------------------------------------')

% Crop results
results.time = results.time(1:save_id-1);

end


%
%
%
%
%
%--------------------------------------------------------------------------
% FUNCTIONS
%--------------------------------------------------------------------------
%
%
%
%
%

%--------------------------------------------------------------------------
% multiple refinements marking all the elements 
%--------------------------------------------------------------------------

function [hmsh, hspace] = uniform_refinement(hmsh, hspace, n_refininements, adaptivity_data)

    if n_refininements >= 1
       
        for i = 1:n_refininements

            if hmsh.nlevels >= adaptivity_data.max_level % check max depth        
                disp('uniform refinement limited by max_level')
                break
            end

            marked = cell(hmsh.nlevels,1); 
            for ilev =1:hmsh.nlevels
                elements = hmsh.active{ilev};  % mark all the active elements
                marked{ilev} = elements;
            end
            [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);

        end
    
    end
end

%--------------------------------------------------------------------------
% adaptivity in space
%--------------------------------------------------------------------------

function [u_n1, udot_n1, hspace, hmsh, est, old_space] = solve_step_adaptive(u_n, udot_n, hspace, hmsh,  ...
                                                dt, a_m, a_f, gamma, lambda, mu, dmu, Cpen, ...
                                                problem_data, adaptivity_data, adaptivity_data_flag, old_space)


    iter = 0;
    while (1)
        
        iter = iter + 1;

        disp(strcat('%%%%%%%%%%%%%%%%% Adaptivity iteration ',num2str(iter),' %%%%%%%%%%%%%%%%%'));
   
        
        %------------------------------------------------------------------
        % solve

        [u_n1, udot_n1, old_space] = generalized_alpha_step(u_n, udot_n, dt, a_m, a_f, gamma, lambda, mu, dmu, Cpen, hspace, hmsh, old_space);
        nel(iter) = hmsh.nel; ndof(iter) = hspace.ndof;


        %------------------------------------------------------------------
        % skip adaptivity
        if adaptivity_data_flag == false
            est = zeros(hmsh.nel, 1);
            disp('No adaptivity')
            break
        end

        %------------------------------------------------------------------
        %estimate
        est = adaptivity_estimate_cahn_hilliard (u_n1, hmsh, hspace, problem_data, adaptivity_data);
        gest(iter) = norm (est);

        %------------------------------------------------------------------
        % stopping criteria

        if (iter == adaptivity_data.num_max_iter)
          disp('Warning: reached the maximum number of iterations')
          solution_data.flag = 2; break
        elseif (hmsh.nlevels > adaptivity_data.max_level)
          disp(strcat('number of levels =',num2str(hmsh.nlevels))) 
          disp('Warning: reached the maximum number of levels')
          solution_data.flag = 3; break             
        end
        
        % stopping criterion
        if hmsh.nlevels == adaptivity_data.max_level
            % check if there are active elements in the coarse levels
            continue_flag = 0;
            for ilev = 1:hmsh.nlevels - 1
                if isempty(hmsh.active{ilev}) == 0
                    continue_flag = 1;
                    break
                end
            end
        end

        if hmsh.nlevels == adaptivity_data.max_level && continue_flag == 0
            disp('Mesh completely refined at the maximum level')
            break
        end

        

        %------------------------------------------------------------------
        % refinement
        % mark
        [marked, num_marked_ref] = adaptivity_mark_cahn_hilliard (est, hmsh, hspace, adaptivity_data);
        % limit the maximum refienment depth
        if hmsh.nlevels == adaptivity_data.max_level
            num_deleted = numel(marked{hmsh.nlevels});
            marked{hmsh.nlevels} = [];
            num_marked_ref = num_marked_ref - num_deleted ;
        end

        % stopping criterion
        if num_marked_ref == 0
            disp('No element is refined')
            break  
        end

        % refine
        [hmsh, hspace, Cref] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
        old_space.modified = 1;
        old_space.space =    [];
        old_space.mesh =     [];
        old_space.mass_mat = [];
        old_space.term3 =    [];
        old_space.term4 =    [];
        old_space.Pen   =    [];
        old_space.pen_rhs =  [];


        % recompute control variables
        [u_n, udot_n] = compute_control_variables_new_mesh(u_n, udot_n, Cref);

   
    end % end loop adaptivity


end


%--------------------------------------------------------------------------
% compute control variables on the new mesh (given the transformation matrix)
%--------------------------------------------------------------------------
   
function [u_n_ref, udot_n_ref] = compute_control_variables_new_mesh(u_n, udot_n, Cref)
       
    u_n_ref = Cref * u_n;       
    udot_n_ref = Cref * udot_n;

end

%--------------------------------------------------------------------------
% coarsening
%--------------------------------------------------------------------------

function [u_n1, udot_n1, hspace, hmsh, old_space] = coarsening_algorithm(est, hmsh, hspace, adaptivity_data,  u_n1, udot_n1, pen_proje, old_space)

    %------------------------------------------------------------------  
    % mark
    [marked, num_marked_coa] = adaptivity_mark_coarsening_cahn_hilliard (est, hmsh, hspace, adaptivity_data);

%     output_file_2 = strcat( 'coarsening_data_debug/Coarsening_data' );
%     save(output_file_2,  'marked', 'est','hmsh', 'hspace')
    
    if num_marked_coa == 0
        old_space.modified = 0;

    else
    %------------------------------------------------------------------
    % coarsening
        hmsh_fine = hmsh;
        hspace_fine = hspace;
        [hmsh, hspace] = adaptivity_coarsen (hmsh, hspace, marked, adaptivity_data);    
        
        if hspace.ndof == hspace_fine.ndof
            old_space.modified = 0;
        else
            % recompute control variables:  (mass+penalty) \ (G * u)
            [u_n1, udot_n1] = compute_control_variables_coarse_mesh(hmsh, hspace, hmsh_fine, hspace_fine, u_n1, udot_n1, pen_proje);

            old_space.modified = 1;
            old_space.space =    [];
            old_space.mesh =     [];
            old_space.mass_mat = [];
            old_space.term3 =    [];
            old_space.term4 =    [];
            old_space.Pen   =    [];
            old_space.pen_rhs =  [];
        end
        clear hmsh_fine hspace_fine

    end
end

%--------------------------------------------------------------------------
% compute control variables on the coarser mesh by means of L2-projection
%--------------------------------------------------------------------------

function [u_n_coa, udot_n_coa] = compute_control_variables_coarse_mesh(hmsh, hspace, hmsh_fine, hspace_fine, u_n, udot_n, pen_proje)

mass_coarse = op_u_v_hier(hspace,hspace,hmsh);

% penalty term (matrix and vector)
bou   = 1:numel(hmsh.mesh_of_level(1).boundaries);
[Pen, ~] = penalty_matrix (hspace, hmsh, bou, pen_proje);

% disp (condest(mass_coarse))

mass_coarse = mass_coarse + Pen;

% disp (condest(mass_coarse))

%G = op_u_v_hier (hspace_fine, hspace_in_finer_mesh(hspace, hmsh, hmsh_fine), hmsh_fine);
%u_n_coa = mass_coarse\(G*u_n);
%udot_n_coa = mass_coarse\(G*udot_n);

rhs_u = op_Gu_hier (hspace_in_finer_mesh(hspace, hmsh, hmsh_fine), hmsh_fine, hspace_fine, u_n);
u_n_coa = mass_coarse\rhs_u;

rhs_udot = op_Gu_hier (hspace_in_finer_mesh(hspace, hmsh, hmsh_fine), hmsh_fine, hspace_fine, udot_n);
udot_n_coa = mass_coarse\rhs_udot;


end

%---------------------------

function rhs = op_Gu_hier (hspace, hmsh_fine, hspace_fine, uhat_fine)



  rhs = zeros (hspace.ndof, 1);

  ndofs = 0;
  for ilev = 1:hmsh_fine.nlevels
    ndofs = ndofs + hspace.ndof_per_level(ilev);
    if (hmsh_fine.nel_per_level(ilev) > 0)
      x = cell (hmsh_fine.rdim, 1);
      for idim = 1:hmsh_fine.rdim
        x{idim} = reshape (hmsh_fine.msh_lev{ilev}.geo_map(idim,:,:), hmsh_fine.msh_lev{ilev}.nqn, hmsh_fine.nel_per_level(ilev));
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

      % 
      b_lev = op_f_v (sp_lev, hmsh_fine.msh_lev{ilev}, u_fine);



      dofs = 1:ndofs;
      rhs(dofs) = rhs(dofs) + hspace.Csub{ilev}.' * b_lev;
    end
  end

end

%--------------------------------------------------------------------------
% One step of generalized alpha method
%--------------------------------------------------------------------------

function [u_n1, udot_n1, old_space] = generalized_alpha_step(u_n, udot_n, dt, a_m, a_f, gamma, ...
    lambda, mu, dmu, Cpen, hspace, hmsh, old_space)

  % Convergence criteria
  n_max_iter = 20;
  tol_rel_res = 1e-10;
  tol_abs_res = 1e-10;

  % Predictor step
  u_n1 = u_n;
  udot_n1 = (gamma-1)/gamma * udot_n; 

  % Newton loop
  for iter = 0:n_max_iter

    % Field at alpha level
    udot_a = udot_n + a_m *(udot_n1-udot_n);
    u_a = u_n + a_f *(u_n1-u_n);

    % Compute the residual (internal)
    [Res_gl, stiff_mat, mass_mat, old_space] = Res_K_cahn_hilliard(hspace, hmsh, lambda, Cpen, u_a, udot_a, mu, dmu, old_space);

    % Convergence check
    if iter == 0
      norm_res_0 = norm(Res_gl);
    end
    norm_res = norm(Res_gl);
    

    if norm_res/norm_res_0 < tol_rel_res
      disp(strcat('iteration n°=',num2str(iter)))
      disp(strcat('norm (abs) residual=',num2str(norm_res)))
      break
    end
    if norm_res<tol_abs_res
      disp(strcat('iteration n°=',num2str(iter)))
      disp(strcat('norm absolute residual=',num2str(norm_res)))
      break
    end
    if iter == n_max_iter
      disp(strcat('Newton reached the maximum number of iterations'))
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
% Canh-Hilliard residual and tangent matrix
%--------------------------------------------------------------------------

function [Res_gl, stiff_mat, mass_mat, old_space] =  Res_K_cahn_hilliard(hspace, hmsh, lambda, Cpen, u_a, udot_a, mu, dmu, old_space)


    % No flux b.c. (essential boundary condition)
    % Set Neumann boundary conditions for non-periodic sides    
    nmnn_bou   = 1:numel(hmsh.mesh_of_level(1).boundaries);


    if old_space.modified == 1
    
      % mass matrix
      mass_mat = op_u_v_hier(hspace,hspace,hmsh);

      % Double well (matrices)
      [term2, term2K] = op_gradfu_gradv_hier_C1(hspace, hmsh, u_a, mu, dmu);   

      % laplacian (matrix)
      lapl_mat = op_laplaceu_laplacev_hier (hspace, hspace, hmsh, lambda);

      % Compute the boundary term (nitsche method)
      bnd_mat = int_boundary_term (hspace, hmsh, lambda, nmnn_bou);
      [pen, pen_rhs] = penalty_matrix (hspace, hmsh, nmnn_bou, Cpen);


      % update old_space
      old_space.modified = 0;
      old_space.space = hspace;
      old_space.mesh = hmsh;
      old_space.mass_mat = mass_mat;
      old_space.lapl_mat = lapl_mat;
      old_space.bnd_mat = bnd_mat;
      old_space.pen   = pen;
      old_space.pen_rhs = pen_rhs;
     
    elseif old_space.modified == 0
      mass_mat =  old_space.mass_mat;
      lapl_mat =  old_space.lapl_mat;
      bnd_mat =  old_space.bnd_mat;
      pen = old_space.pen;
      pen_rhs =  old_space.pen_rhs;
      [term2, term2K] = op_gradfu_gradv_hier_C1(old_space.space, old_space.mesh, u_a, mu, dmu);

    end

    %----------------------------------------------------------------------
    % Residual
    Res_gl = mass_mat*udot_a   + lapl_mat*u_a + term2*u_a;

    % Tangent stiffness matrix (mass is not considered here)
    stiff_mat =  lapl_mat + term2 + term2K ;

    % nitsche method to impose Neumann B.C.
    Res_gl = Res_gl - (bnd_mat + bnd_mat.') * u_a + pen*u_a - pen_rhs;
    stiff_mat = stiff_mat - (bnd_mat + bnd_mat.') + pen;


    
end

%--------------------------------------------------------------------------
% Boundary term, \int_\Gamma (\Delta u) (\partial v / \partial n)
%--------------------------------------------------------------------------

function [A] = int_boundary_term (hspace, hmsh,  lambda, nmnn_sides)

  if (~isempty(nmnn_sides))

      A = spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);

    
      boundaries = hmsh.mesh_of_level(1).boundaries;
      Nbnd = cumsum ([0, boundaries.nsides]);
      last_dof = cumsum (hspace.ndof_per_level);
        
      for iref = nmnn_sides
        iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
        for ilev = 1:hmsh.boundary.nlevels
          patch_bnd_shifting = cumsum ([0 hmsh.boundary.mesh_of_level(ilev).nel_per_patch]);
          patch_shifting = cumsum ([0 hmsh.mesh_of_level(ilev).nel_per_patch]);
          dofs_on_lev = 1:last_dof(ilev);
    
          for ii = 1:numel(iref_patch_list)
            iptc_bnd = iref_patch_list(ii);
            iptc = boundaries(iref).patches(ii);
            iside = boundaries(iref).faces(ii);
            elems_patch = patch_bnd_shifting(iptc_bnd)+1:patch_bnd_shifting(iptc_bnd+1);
            [~, ~, elements] = intersect (hmsh.boundary.active{ilev}, elems_patch);
    
            if (~isempty (elements))
              msh_patch = hmsh.mesh_of_level(ilev).msh_patch{iptc};
    
              msh_side = msh_eval_boundary_side (msh_patch, iside, elements);
              msh_side_from_interior = msh_boundary_side_from_interior (msh_patch, iside);
    
              sp_bnd = hspace.space_of_level(ilev).sp_patch{iptc}.constructor (msh_side_from_interior);
              msh_side_from_interior_struct = msh_evaluate_element_list (msh_side_from_interior, elements);
              sp_bnd = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', false, 'gradient', true, 'laplacian', true);
    
              for idim = 1:hmsh.rdim
                x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
              end
              coe_side = lambda (x{:});

              [Cpatch, Cpatch_cols_lev] = sp_compute_Cpatch (hspace.space_of_level(ilev), iptc);
              [~,Csub_rows,Cpatch_cols] = intersect (hspace.Csub_row_indices{ilev}, Cpatch_cols_lev);
              Caux = Cpatch(:,Cpatch_cols) * hspace.Csub{ilev}(Csub_rows,:);


              tmp = op_gradv_n_laplaceu(sp_bnd ,sp_bnd ,msh_side, coe_side);   

              A(dofs_on_lev,dofs_on_lev) = A(dofs_on_lev, dofs_on_lev) + Caux.' * tmp * Caux;
             end
           end
        end
      end

  else
    A = [];
  end
end

%--------------------------------------------------------------------------
% penalty term
%--------------------------------------------------------------------------

function [A, rhs] = penalty_matrix (hspace, hmsh, nmnn_sides, Cpen)


  A = spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);
  rhs = zeros (hspace.ndof, 1);

  boundaries = hmsh.mesh_of_level(1).boundaries;
  Nbnd = cumsum ([0, boundaries.nsides]);
  last_dof = cumsum (hspace.ndof_per_level);
    
  for iref = nmnn_sides
    iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    for ilev = 1:hmsh.boundary.nlevels
      patch_bnd_shifting = cumsum ([0 hmsh.boundary.mesh_of_level(ilev).nel_per_patch]);
      patch_shifting = cumsum ([0 hmsh.mesh_of_level(ilev).nel_per_patch]);
      dofs_on_lev = 1:last_dof(ilev);

      for ii = 1:numel(iref_patch_list)
        iptc_bnd = iref_patch_list(ii);
        iptc = boundaries(iref).patches(ii);
        iside = boundaries(iref).faces(ii);
        elems_patch = patch_bnd_shifting(iptc_bnd)+1:patch_bnd_shifting(iptc_bnd+1);
        [~, ~, elements] = intersect (hmsh.boundary.active{ilev}, elems_patch);

        if (~isempty (elements))
          msh_patch = hmsh.mesh_of_level(ilev).msh_patch{iptc};

          msh_side = msh_eval_boundary_side (msh_patch, iside, elements);
          msh_side_from_interior = msh_boundary_side_from_interior (msh_patch, iside);

          sp_bnd = hspace.space_of_level(ilev).sp_patch{iptc}.constructor (msh_side_from_interior);
          msh_side_from_interior_struct = msh_evaluate_element_list (msh_side_from_interior, elements);
          sp_bnd = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', true, 'gradient', true);

%           for idim = 1:hmsh.rdim
%             x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
%           end
%           
          % For simplicity I replaced charlen by 2^(-level)
%           coeff_at_qnodes =  coeff_at_qnodes * Cpen ./ msh_side.charlen;
          coeff_at_qnodes = Cpen * ones(msh_side.nqn, msh_side.nel)/ 2^(-ilev);

          [Cpatch, Cpatch_cols_lev] = sp_compute_Cpatch (hspace.space_of_level(ilev), iptc);
          [~,Csub_rows,Cpatch_cols] = intersect (hspace.Csub_row_indices{ilev}, Cpatch_cols_lev);
          Caux = Cpatch(:,Cpatch_cols) * hspace.Csub{ilev}(Csub_rows,:);
          
          
          tmp = op_gradu_n_gradv_n(sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
          tmp = Caux.' * tmp * Caux;

          A(dofs_on_lev,dofs_on_lev) = A(dofs_on_lev, dofs_on_lev) + tmp;
%           rhs(dofs_on_lev) = rhs(dofs_on_lev) + Caux.' * (tmp_rhs);
        end
      end
    end
  end



end

%--------------------------------------------------------------------------
% double-well function
%--------------------------------------------------------------------------

function [A, B] = op_gradfu_gradv_hier_C1(hspace, hmsh, uhat, f, df)

  A = spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);
  B = spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);
  patch_list = 1:hmsh.npatch;

  ndofs_u = 0;
  ndofs_v = 0;
  last_dof = cumsum (hspace.ndof_per_level);



  for ilev = 1:hmsh.nlevels
    ndofs_u = ndofs_u + hspace.ndof_per_level(ilev);
    ndofs_v = ndofs_v + hspace.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)
      msh_lev = msh_restrict_to_patches (hmsh.msh_lev{ilev}, patch_list);

      if (msh_lev.nel > 0)
        x = cell(msh_lev.rdim,1);
        for idim = 1:hmsh.rdim
          x{idim} = reshape (msh_lev.geo_map(idim,:,:), msh_lev.nqn, msh_lev.nel);
        end
        spu_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), msh_lev, 'value', true, 'gradient', true);
        spu_lev = change_connectivity_localized_Csub (spu_lev, hspace, ilev);

        dofs_u = 1:ndofs_u;
        dofs_v = 1:ndofs_v;


        % Evaluate the field and its gradient at the Gaussian points
        uhat_lev = hspace.Csub{ilev}*uhat(1:last_dof(ilev));
        utemp = sp_eval_msh (uhat_lev, spu_lev, msh_lev, {'value', 'gradient'});
        u = utemp{1};
        gradu = utemp{2};
        

        % double-well
        coeffs_A = f(u);
        A_lev = op_gradu_gradv (spu_lev, spu_lev, msh_lev, coeffs_A);
        A(dofs_v,dofs_u) = A(dofs_v,dofs_u) + hspace.Csub{ilev}.' * A_lev * hspace.Csub{ilev};

        coeffs_B = df(u);
        coeffs_Bv = gradu;
        for idim = 1:hmsh.ndim
          coeffs_Bv(idim,:,:) = coeffs_Bv(idim,:,:) .* reshape(coeffs_B, 1, size(coeffs_B,1), size(coeffs_B,2));
        end
        B_lev = op_vel_dot_gradu_v (spu_lev, spu_lev, msh_lev, coeffs_Bv).';     
        B(dofs_v,dofs_u) = B(dofs_v,dofs_u) + hspace.Csub{ilev}.' * B_lev * hspace.Csub{ilev};


      end
    end
  end

end

%--------------------------------------------------------------------------
% save and plot results
%--------------------------------------------------------------------------

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
