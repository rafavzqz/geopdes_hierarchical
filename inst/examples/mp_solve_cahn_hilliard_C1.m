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

function [geometry, msh, space, results] = mp_solve_cahn_hilliard_C1 (problem_data, method_data, save_info)

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

%%-------------------------------------------------------------------------
% Parameters for the double-well function (to be given in problem_data)
alpha = 1; 
beta = 1;
mu = @(x) 3 * alpha * x.^2 - beta;
dmu = @(x) 6 * alpha * x;

%%-------------------------------------------------------------------------
% Construct geometry structure, and information for interfaces and boundaries
[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

msh = cell (1, npatch); 
sp = cell (1, npatch);
for iptc = 1:npatch
% Define the refined mesh, with tensor product structure
  [knots{iptc}, zeta{iptc}] = ...
         kntrefine (geometry(iptc).nurbs.knots, nsub-1, degree, regularity);

% Compute the quadrature rule
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (zeta{iptc}, rule);
  msh{iptc} = msh_cartesian (zeta{iptc}, qn, qw, geometry(iptc));

% Evaluate the discrete space basis functions in the quadrature points
  sp{iptc} = sp_bspline (knots{iptc}, degree, msh{iptc});
end

[edges, vertices] = vertices_struct (geometry, interfaces, boundaries, boundary_interfaces);
msh = msh_multipatch (msh, boundaries);
space = sp_multipatch_C1 (sp, msh, geometry, edges, vertices);
clear sp

%%-------------------------------------------------------------------------
% Generalized alpha parameters
a_m = .5*((3-rho_inf_gen_alpha)/(1+rho_inf_gen_alpha));
a_f =  1/(1+rho_inf_gen_alpha);
gamma =  .5 + a_m - a_f;

%%-------------------------------------------------------------------------
% No flux b.c. (essential boundary condition)
% Set Neumann boundary conditions for non-periodic sides
nmnn_bou   = 1:numel(boundaries);

%%-------------------------------------------------------------------------
% Precompute matrices

% Compute the mass matrix
mass_mat = op_u_v_mp (space,space,msh);

% Compute the laplace matrix
lapl_mat = op_laplaceu_laplacev_mp (space, space, msh, lambda);

% Compute the boundary term
bnd_mat = int_boundary_term (space, msh, lambda, nmnn_bou);

% Compute the penalty matrix
[Pen, pen_rhs] = penalty_matrix (space, msh, nmnn_bou, Cpen_nitsche);


%%-------------------------------------------------------------------------
% Initial conditions
time = 0;
if (exist('fun_u', 'var') && ~isempty(fun_u))
  rhs = op_f_v_mp (space, msh, fun_u);
  u_n = (mass_mat + Cpen_projection/Cpen_nitsche * Pen)\rhs;
else
  u_n = zeros(space.ndof, 1);
end

if (exist('fun_udot', 'var') && ~isempty(fun_udot))
  rhs = op_f_v_mp (space, msh, fun_udot);
  udot_n = (mass_mat + Cpen_projection/Cpen_nitsche * Pen)\rhs;
else
  udot_n = zeros(space.ndof, 1);
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
                    mass_mat, lapl_mat, bnd_mat, Pen, pen_rhs, space, msh);



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
    if (time >= save_info(save_id))
      save_id = save_id + 1;
      results.u(:,save_id) = u_n;
      results.udot(:,save_id) = udot_n;
      results.time(save_id) = time;
      if (save_id > length(save_info))
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
  [term2, term2K] = op_gradfu_gradv_mp (space, msh, u_a, mu, dmu);    
 
  % Residual
  Res_gl = mass_mat*udot_a + term2*u_a  + lapl_mat*u_a;

  % Tangent stiffness matrix (mass is not considered here)
  stiff_mat = term2 + term2K + lapl_mat;

  % in case of neumann BC, add boundary terms
  if (~isempty(bnd_mat))
    Res_gl = Res_gl - (bnd_mat + bnd_mat.') * u_a + Pen*u_a - pen_rhs;
    stiff_mat = stiff_mat - (bnd_mat + bnd_mat.') + Pen;
  end
end

%--------------------------------------------------------------------------
% Boundary term, \int_\Gamma (\Delta u) (\partial v / \partial n)
%--------------------------------------------------------------------------

function [A] = int_boundary_term (space, msh,  lambda, nmnn_sides)

  if (~isempty(nmnn_sides))

    A =  spalloc (space.ndof, space.ndof, 3*space.ndof);

    for iref = nmnn_sides
      for bnd_side = 1:msh.boundaries(iref).nsides  
        iptc = msh.boundaries(iref).patches(bnd_side);
        iside = msh.boundaries(iref).faces(bnd_side);

        msh_side = msh_eval_boundary_side (msh.msh_patch{iptc}, iside ) ;
        msh_side_int = msh_boundary_side_from_interior (msh.msh_patch{iptc}, iside ) ;

        sp_side = space.sp_patch{iptc}.constructor (msh_side_int);
        sp_side = sp_precompute ( sp_side , msh_side_int , 'gradient' , true, 'laplacian', true );

        for idim = 1:msh.rdim
          x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
        end
        coe_side = lambda (x{:});
        tmp = op_gradv_n_laplaceu(sp_side ,sp_side ,msh_side, coe_side);

        [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, iptc);
        A(Cpatch_cols,Cpatch_cols) = ...
              A(Cpatch_cols,Cpatch_cols) + Cpatch.' * tmp * Cpatch;
      end
    end

  else
    A = [];
  end
end

%--------------------------------------------------------------------------
% penalty term
%--------------------------------------------------------------------------

function [P, rhs] = penalty_matrix (space, msh, nmnn_sides, pen)

  P =  spalloc (space.ndof, space.ndof, 3*space.ndof);
  rhs = zeros(space.ndof,1);

  for iref = nmnn_sides
    for bnd_side = 1:msh.boundaries(iref).nsides
      iptc = msh.boundaries(iref).patches(bnd_side);
      iside = msh.boundaries(iref).faces(bnd_side);

      msh_side = msh_eval_boundary_side (msh.msh_patch{iptc}, iside);
      msh_side_from_interior = msh_boundary_side_from_interior (msh.msh_patch{iptc}, iside);

      sp_bnd = space.sp_patch{iptc}.constructor (msh_side_from_interior);
      sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', false, 'gradient', true);

      coe_side = pen .* msh_side.charlen; 
      tmp = op_gradu_n_gradv_n(sp_bnd, sp_bnd, msh_side, coe_side);
      
      [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, iptc);
      P(Cpatch_cols,Cpatch_cols) = P(Cpatch_cols,Cpatch_cols) + ...
        Cpatch.' * (tmp) * Cpatch;

      %rhs = rhs + rhs_pen;
    end
  end

end


% function [mass_pen, rhs_pen] = penalty_grad (space, msh, i_side, pen)
% 
% msh_side = msh_eval_boundary_side (msh, i_side ) ;
% msh_side_int = msh_boundary_side_from_interior (msh, i_side ) ;
% sp_side = space.constructor ( msh_side_int) ;
% sp_side = sp_precompute ( sp_side , msh_side_int , 'gradient', true );
% 
% 
% coe_side = pen .* msh_side.charlen; 
% mass_pen = op_gradu_n_gradv_n(sp_side, sp_side, msh_side, coe_side);
% 
% 
% rhs_pen  = zeros(space.ndof,1); % no flux
% 
% end

%--------------------------------------------------------------------------
% check flux through the boundaries
%--------------------------------------------------------------------------

% function flux = check_flux_phase_field(space, msh, uhat, uhat0)
% 
% sides = [1,2,3,4]; % all the boundaries
% flux = 0;
% 
% for iside=1:length(sides)   
% 
%     msh_side = msh_eval_boundary_side (msh, sides(iside) ) ;
%     msh_side_int = msh_boundary_side_from_interior (msh, sides(iside) ) ;
%     sp_side = space.constructor ( msh_side_int) ;
%     sp_side = sp_precompute ( sp_side , msh_side_int , 'gradient' , true );
% 
% 
%     gradu = sp_eval_msh (uhat-uhat0, sp_side, msh_side, 'gradient');
%    
% 
%     
%     valu = zeros(sp_side.ncomp, size(msh_side.quad_weights,1), size(msh_side.quad_weights,2));
%     for idim = 1:msh.rdim
%         valu = valu + (gradu(idim,:,:) .* msh_side.normal(idim,:,:));
%     end
% 
% 
%     w =msh_side.quad_weights .* msh_side.jacdet;
%     err_elem = sum (reshape (valu, [msh_side.nqn, msh_side.nel]) .* w);
%     err  = sum (err_elem);
% 
% 
%     flux = flux + err;
% 
% end
% 
% 
% 
% end

