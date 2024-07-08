% PHYSICAL DATA OF THE PROBLEM
clear problem_data

% PHYSICAL DATA OF THE PROBLEM
problem_data.geo_name = 'geo_open_quasisphere_5p_ASG1.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];
problem_data.weak_drchlt_sides = [];

% Physical parameters
problem_data.c_diff  = @(x, y, z) ones(size(x));

% Source and boundary terms
problem_data.f = @(x, y, z) ones(size(x));
problem_data.g = @(x, y, z, ind) zeros(size(x));
problem_data.h = @(x, y, z, ind) zeros(size(x));

% DISCRETIZATION PARAMETERS
clear method_data
deg = 4;
method_data.degree      = [deg deg];     % Degree of the splines
method_data.regularity  = [deg-2 deg-2]; % Regularity of the splines
method_data.nsub_coarse = [4 4];         % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];         % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nquad       = [deg+1 deg+1]; % Points for the Gaussian quadrature rule
method_data.space_type  = 'standard';
method_data.truncated   = 1;            % 0: False, 1: True
method_data.interface_regularity = 1;

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag          = 'elements';
adaptivity_data.C0_est        = 1.0;
adaptivity_data.mark_param    = .25;
adaptivity_data.mark_strategy = 'GERS';
adaptivity_data.max_level     = 10;
adaptivity_data.max_ndof      = 5000;
adaptivity_data.num_max_iter  = 5;
adaptivity_data.max_nel       = 5000;
adaptivity_data.tol           = 1e-5;

% GRAPHICS
plot_data.print_info = true;
plot_data.plot_hmesh = false;
plot_data.plot_discrete_sol = false;

[geometry, hmsh, hspace, u, solution_data] = adaptivity_bilaplace_mp_C1 (problem_data, method_data, adaptivity_data, plot_data);
npts = [15 15];
sp_plot_solution (u, hspace, geometry, npts); shading interp; axis equal
