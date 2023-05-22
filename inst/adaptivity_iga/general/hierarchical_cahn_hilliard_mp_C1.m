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
for iopt  = 1:numel (initial_conditions)
  eval ([data_names{iopt} '= initial_conditions.(data_names{iopt});']);
end


%%-------------------------------------------------------------------------
% Parameters for the double-well function (to be given in problem_data)
alpha = 1; 
beta = 1;
mu = @(x) 3 * alpha * x.^2 - beta;
dmu = @(x) 6 * alpha * x;

%%-------------------------------------------------------------------------
% Construct geometry structure, and information for interfaces and boundaries

% Initialization of the most coarse level of the hierarchical mesh and space
[hmsh, hspace, geometry] = adaptivity_initialize_laplace_mp_C1 (problem_data, method_data);

% Refine the mesh up to a predefined level
n_refininements = adaptivity_data.max_level-1; % number of uniform refinements
[hmsh, hspace] = uniform_refinement(hmsh, hspace, n_refininements, adaptivity_data);

%%-------------------------------------------------------------------------
% Generalized alpha parameters
a_m = .5*((3-rho_inf_gen_alpha)/(1+rho_inf_gen_alpha));
a_f =  1/(1+rho_inf_gen_alpha);
gamma =  .5 + a_m - a_f;

%%-------------------------------------------------------------------------
% No flux b.c. (essential boundary condition)
% Set Neumann boundary conditions for non-periodic sides

nmnn_bou   = 1:numel(hmsh.mesh_of_level(1).boundaries);

%%-------------------------------------------------------------------------
% Precompute matrices

% Compute the mass matrix
mass_mat = op_u_v_hier (hspace,hspace,hmsh);

% Compute the laplace matrix
lapl_mat = op_laplaceu_laplacev_hier (hspace, hspace, hmsh, lambda);

% Compute the boundary term
bnd_mat = int_boundary_term (hspace, hmsh, lambda, nmnn_bou);

% Compute the penalty matrix
[Pen, pen_rhs] = penalty_matrix (hspace, hmsh, nmnn_bou, Cpen_nitsche);


%%-------------------------------------------------------------------------
% Initial conditions
time = 0;
if (exist('fun_u', 'var') && ~isempty(fun_u))
  rhs = op_f_v_hier (hspace, hmsh, fun_u);
  u_n = (mass_mat + Cpen_projection/Cpen_nitsche * Pen)\rhs;
else
  u_n = zeros(hspace.ndof, 1);
end

if (exist('fun_udot', 'var') && ~isempty(fun_udot))
  rhs = op_f_v_hier (space, msh, fun_udot);
  udot_n = (mass_mat + Cpen_projection/Cpen_nitsche * Pen)\rhs;
else
  udot_n = zeros(hspace.ndof, 1);
end



%%-------------------------------------------------------------------------
% Initialize structure to store the results
save_id = 1;
results.u = zeros(length(u_n), length(save_info)+1);
results.udot = zeros(length(u_n), length(save_info)+1);
results.time = zeros(length(save_info)+1,1);
flag_stop_save = false;

% Save initial conditions
results.u(:,1) = u_n;
results.udot(:,1) = udot_n;
results.time(1) = time;

%%-------------------------------------------------------------------------
% Loop over time steps
while time < Time_max
  disp('----------------------------------------------------------------')
  disp(strcat('time step t=',num2str(time)))

  [u_n1, udot_n1] = generalized_alpha_step(u_n, udot_n, dt, a_m, a_f, gamma, mu, dmu, ...
                    mass_mat, lapl_mat, bnd_mat, Pen, pen_rhs, hspace, hmsh);


  % Time step update
  time = time + dt;
  u_n = u_n1;
  udot_n = udot_n1;

  % Check max time
  if (time + dt > Time_max)
    dt = Time_max-time;
  end

  % Store results
  if (flag_stop_save == false)
    if (time >= save_info.time_save(save_id))
      save_id = save_id + 1;
      results.u(:,save_id) = u_n;
      results.udot(:,save_id) = udot_n;
      results.time(save_id) = time;
      if (save_id > length(save_info.time_save))
        flag_stop_save = true;
      end
    end
  end
end

disp('----------------------------------------------------------------')
disp(strcat('END ANALYSIS t=',num2str(time)))
disp('----------------------------------------------------------------')

% Crop results
results.u = results.u(:,1:save_id);
results.udot = results.udot(:,1:save_id);
results.time = results.time(1:save_id);

end

%--------------------------------------------------------------------------
% FUNCTIONS
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% One step of generalized alpha-method
%--------------------------------------------------------------------------

function [u_n1, udot_n1] = generalized_alpha_step(u_n, udot_n, dt, a_m, a_f, gamma, mu, dmu, ...
                       mass_mat, lapl_mat, bnd_mat, Pen, pen_rhs, space, msh)

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
    [Res_gl, stiff_mat] = Res_K_cahn_hilliard (space, msh, ...
                          mass_mat, lapl_mat, bnd_mat, Pen, pen_rhs, ...
                          u_a, udot_a, mu, dmu);

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

function [Res_gl, stiff_mat] = Res_K_cahn_hilliard(space, msh, ...
                          mass_mat, lapl_mat, bnd_mat, Pen, pen_rhs, u_a, udot_a, mu, dmu)

  % Double well (matrices)
  [term2, term2K] = op_gradfu_gradv_hier_C1(space, msh, u_a, mu, dmu);    
 
  % Residual
  Res_gl = mass_mat*udot_a   + lapl_mat*u_a + term2*u_a;

  % Tangent stiffness matrix (mass is not considered here)
  stiff_mat =  lapl_mat + term2 + term2K ;

  % in case of neumann BC, add boundary terms
  if (~isempty(bnd_mat))
    Res_gl = Res_gl - (bnd_mat + bnd_mat.') * u_a + Pen*u_a - pen_rhs;
    stiff_mat = stiff_mat - (bnd_mat + bnd_mat.') + Pen;
  end
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
        utemp = sp_eval_msh (uhat(dofs_u), spu_lev, msh_lev, {'value', 'gradient'});
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





function [A, B] = op_gradfu_gradv (space, msh, uhat, f, df)

  for idim = 1:msh.ndim
    size1 = size (space.sp_univ(idim).connectivity);
    if (size1(2) ~= msh.nel_dir(idim))
      error ('The discrete space is not associated to the mesh')
    end
  end

  A = spalloc (space.ndof, space.ndof, 3*space.ndof);
  B = spalloc (space.ndof, space.ndof, 6*space.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'gradient', true);

    % Evaluate the field and its gradient at the Gaussian points
    utemp = sp_eval_msh (uhat, sp_col, msh_col, {'value', 'gradient'});
    u = utemp{1};
    gradu = utemp{2};

    % Polynomial formulation for the double-well
    coeffs_A = f(u);
    A = A + op_gradu_gradv (sp_col, sp_col, msh_col, coeffs_A);

    coeffs_B = df(u);
    coeffs_Bv = gradu;
    for idim = 1:msh.ndim
      coeffs_Bv(idim,:,:) = coeffs_Bv(idim,:,:) .* reshape(coeffs_B, 1, size(coeffs_B,1), size(coeffs_B,2));
    end
    B = B + op_vel_dot_gradu_v (sp_col, sp_col, msh_col, coeffs_Bv).';
  end
end

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


