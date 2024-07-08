% PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_3patch_ASG1.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];
problem_data.weak_drchlt_sides = [1 2 3 4 5 6];

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));

% Source and boundary terms
C = 20; 
P = [17/3+0.15, 2+0.075];
alpha = 4/6;
problem_data.f = @(x,y) -2*alpha*((x-P(1)).^2+(y-P(2)).^2).^(alpha-2).*...
                       (((x-P(1)).^2+(y-P(2)).^2)+2*(alpha-1)*(x-P(1)).^2)...
                        -2*alpha*((x-P(1)).^2+(y-P(2)).^2).^(alpha-2).*...
                       (((x-P(1)).^2+(y-P(2)).^2)+2*(alpha-1)*(y-P(2)).^2);
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) ((x-P(1)).^2+(y-P(2)).^2).^alpha;

% Exact solution (optional)
problem_data.uex =@(x,y) ((x-P(1)).^2+(y-P(2)).^2).^alpha;
problem_data.graduex = @(x,y) 2*alpha*cat (1, ...
            reshape (((x-P(1)).^2+(y-P(2)).^2).^(alpha-1).*(x-P(1)), [1, size(x)]), ...
            reshape (((x-P(1)).^2+(y-P(2)).^2).^(alpha-1).*(y-P(2)), [1, size(x)]));

% DISCRETIZATION PARAMETERS
clear method_data
deg = 4;
method_data.degree      = [deg deg];     % Degree of the splines
method_data.regularity  = [deg-2 deg-2]; % Regularity of the splines
method_data.nsub_coarse = [2 2];         % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];         % Number of subdivisions for each refinement
method_data.nquad       = [deg+1 deg+1]; % Points for the Gaussian quadrature rule
method_data.Cpen = 10 * (min(method_data.degree) + 1);
method_data.space_type  = 'simplified';  % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 1;             % 0: False, 1: True
method_data.interface_regularity = 1;

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag          = 'elements';
adaptivity_data.C0_est        = 1.0;
adaptivity_data.mark_param    = .2;
adaptivity_data.mark_strategy = 'GERS';
adaptivity_data.max_level     = 10;
adaptivity_data.max_ndof      = 5000;
adaptivity_data.num_max_iter  = 6;
adaptivity_data.max_nel       = 5000;
adaptivity_data.tol           = 1e-5;

% GRAPHICS
plot_data.print_info = true;
plot_data.plot_hmesh = false;
plot_data.plot_discrete_sol = false;

[geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace_mp_C1 (problem_data, method_data, adaptivity_data, plot_data);

% EXPORT VTK FILE
npts = [51 51];
output_file = 'laplace_adaptivity_3patches.pvd';
sp_to_vtk (u, hspace, geometry, npts, output_file, {'solution', 'gradient', 'laplacian'}, {'value', 'gradient', 'laplacian'})

sp_plot_solution (u, hspace, geometry, npts); shading interp
