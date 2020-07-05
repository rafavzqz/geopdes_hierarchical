% EX_KIRCHHOFF_RECTANGULAR_PLATE_CLAMPED_ADAPTIVITY: solve bilaplacian in a clamped rectangular plate.
% PHYSICAL DATA OF THE PROBLEM
clear classes
clear problem_data
% Geometry definition (unrefined geometry)
base = 1.0; height = 1.0;
p11 =[0 0]; p12 =[base 0]; p21 =[0 height]; p22 =[base height];
srf = nrb4surf (p11,p12,p21,p22);

problem_data.geo_name = srf;

% Elasticity constants
E  =  12.0;                         % Young modulus   
nu = 0.0;                         % Poisson modulus
thickness = 1.0;                         % Thickness of the plate
D  = E*thickness^3.0/(12.0*(1.0-nu*nu));   % Flexural rigidity of the plate
p = -1.0;                           % Distributed load
max_iteration = 400;   % number of iteration for computing the exact solution of problem (infinite double sum of sin)

% Boundary conditions
problem_data.simply_supported_sides = [1 2 3 4];
problem_data.clamped_sides = [];
problem_data.nmnn_sides = [];
problem_data.press_sides = [];
problem_data.drchlt_sides = [];
problem_data.prescribed_moment_sides = [];
problem_data.h = @(x, y, ind)  zeros(size(x));

% Physical parameters
problem_data.D = @(x, y) D*ones(size(x));
% Source term  
problem_data.f = @(x, y) p*ones(size(x));
% CHOICE OF THE DISCRETIZATION PARAMETERS
%[grad_x, grad_y] = rectangular_plate_uniform_load_analytical_gradient (x, y, a, b, D, q0, max_iteration, 'x')
problem_data.uex     = @(x, y) rectangular_plate_uniform_load_analytical_solution (x, y, base, height, D, p, max_iteration);

problem_data.graduex = @(x, y) cat (1, ...
                       reshape (rectangular_plate_uniform_load_analytical_gradient (x, y, base, height, D, p, max_iteration, 'x'), [1, size(x)]), ...
                       reshape (rectangular_plate_uniform_load_analytical_gradient (x, y, base, height, D, p, max_iteration, 'y'), [1, size(x)]));

problem_data.hessuex = @(x, y) cat (1, ...
                       reshape (rectangular_plate_uniform_load_analytical_hessian (x, y, base, height, D, p, max_iteration, 'xx'), [1, size(x)]), ...
                       reshape (rectangular_plate_uniform_load_analytical_hessian (x, y, base, height, D, p, max_iteration, 'xy'), [1, size(x)]), ...
                       reshape (rectangular_plate_uniform_load_analytical_hessian (x, y, base, height, D, p, max_iteration, 'yx'), [1, size(x)]), ...
                       reshape (rectangular_plate_uniform_load_analytical_hessian (x, y, base, height, D, p, max_iteration, 'yy'), [1, size(x)]));

n = 2;
deg = 4;

clear method_data
method_data.degree      = [deg deg];        % Degree of the splines
method_data.regularity  = [deg-1 deg-1];        % Regularity of the splines
method_data.nsub_coarse = [n n];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];        % Number of subdivisions for each refinement
method_data.nquad       = [deg+1 deg+1];        % Points for the Gaussian quadrature rule
method_data.space_type  = 'simplified'; % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 0;            % 0: False, 1: True

clear adaptivity_data
adaptivity_data.mark_param = .5;
adaptivity_data.mark_strategy = 'MS';
% adaptivity_data.mark_strategy = 'GR';
adaptivity_data.flag = 'elements';
% adaptivity_data.flag = 'functions';
adaptivity_data.num_max_iter = 5;
adaptivity_data.max_nel = 5000;
adaptivity_data.tol = 1e-10;

geometry = geo_load (problem_data.geo_name);
[knots, zeta] = kntrefine (geometry.nurbs.knots, method_data.nsub_coarse-1, method_data.degree, method_data.regularity);

rule     = msh_gauss_nodes (method_data.nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh   = msh_cartesian (zeta, qn, qw, geometry, 'der2', true);
space = sp_bspline (knots, method_data.degree, msh);
hmsh     = hierarchical_mesh (msh, method_data.nsub_refine);
hspace   = hierarchical_space (hmsh, space, method_data.space_type, method_data.truncated);

nel = zeros (1, adaptivity_data.num_max_iter+1); ndof = nel; gest = nel+NaN;
dof_vector = zeros(1, adaptivity_data.num_max_iter+1);
if (isfield (problem_data, 'hessuex'))
    err_h2 = gest;
end

%% ADAPTIVE LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter = 0;
while (iter < adaptivity_data.num_max_iter)
  iter = iter + 1;
  
  ndof = hspace.ndof;
  
  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',iter);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% SOLVE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('SOLVE:')
  % Assemblying stiffness matrix
  % CALL TO THE SOLVER
  u = adaptivity_solve_bilaplace_gradgrad_2d_iso(hspace, hmsh, problem_data);
  
% POST-PROCESSING
% EXPORT TO PARAVIEW
% output_file = sprintf('Kirchhoff_Rectangular_Plate_Simply_Supported_Adaptivity %d',iter );
% vtk_pts = {linspace(0, 1, 31), linspace(0, 1, 31)};
% % fprintf ('The result is saved in the file %s \n \n', output_file);
% % sp_to_vtk (u, hspace, geometry, vtk_pts, output_file, 'u')
% 
% eu = sp_eval(u, hspace, geometry, vtk_pts);
% max_displacement = min (eu(:));
% fprintf ('Computed solution, max. displacement = %e \n', max_displacement);
% analytical_solution = ...
%     rectangular_plate_analytical_solution (base, height, thickness, E, nu, p);
% fprintf ('Analytical solution, max. displacement = %e \n', analytical_solution);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% ESTIMATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  estimator = adaptivity_estimate_bilaplace (u, hmsh, hspace, problem_data, adaptivity_data);

  if (isfield (problem_data, 'hessuex'))
%     [err_h1(iter), err_l2(iter), err_h1s(iter)] = sp_h1_error (hspace, hmsh, u, problem_data.uex, problem_data.graduex);
    [err_h2(iter), ~, ~, ~, ~, ~, ~, ~, ~, ~] = sp_h2_error (hspace, hmsh, u, ...
                                            problem_data.uex, problem_data.graduex, problem_data.hessuex);
    dof_vector(iter) = ndof;
  end

%   space_bubble = space_bubble_function (hmsh, 'plate');
%   estimator = compute_bubble_estimator_bilaplacian (hspace, space_bubble, hmsh, problem_data.f, u, problem_data.D);
  
  gest(iter) = norm (estimator);
%   gest = zeros (1, adaptivity_data.num_max_iter)+NaN;  
%   if (gest(iter) < adaptivity_data.tol)
%     disp('Success: The error estimation reached the desired tolerance'); 
%     solution_data.flag = 1; break
%   elseif (iter == adaptivity_data.num_max_iter)
%     disp('Warning: reached the maximum number of iterations')
%     solution_data.flag = 2; break
%   end
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% MARK
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('MARK:')
  [marked, num_marked] = adaptivity_mark (estimator, hmsh, hspace, adaptivity_data);
  fprintf('%d %s marked for refinement \n', num_marked, adaptivity_data.flag);
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% REFINE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('REFINE:')
    [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
  fprintf('\n');
  
end

u = adaptivity_solve_bilaplace_gradgrad_2d_iso(hspace, hmsh, problem_data);
% output_file = sprintf('Kirchhoff_Rectangular_Plate_Simply_Supported_Adaptivity %d',iter+1 );
% vtk_pts = {linspace(0, 1, 31), linspace(0, 1, 31)};
% fprintf ('The result is saved in the file %s \n \n', output_file);
% sp_to_vtk (u, hspace, geometry, vtk_pts, output_file, 'u')
% hmsh_plot_cells(hmsh);

estimator = adaptivity_estimate_bilaplace (u, hmsh, hspace, problem_data, adaptivity_data);

% space_bubble = space_bubble_function (hmsh, 'plate');
% estimator = compute_bubble_estimator_bilaplacian (hspace, space_bubble, hmsh, problem_data.f, u, problem_data.D);

gest(iter+1) = norm (estimator);

if (isfield (problem_data, 'hessuex'))
    [err_h2(iter+1), ~, ~, ~, ~, ~, ~, ~, ~, ~] = sp_h2_error (hspace, hmsh, u, ...
                                        problem_data.uex, problem_data.graduex, problem_data.hessuex);
    dof_vector(iter+1) = hspace.ndof;
end

disp('END:')

