% PHYSICAL DATA OF THE PROBLEM
clear problem_data

problem_data.geo_name = 'geo_paraboloid_ASG1.txt';

problem_data.drchlt_sides = [1 2 3 4 5 6 7 8];
problem_data.drchlt_components = {[1 2 3] [1 2 3] [1 2 3] [1 2 3] [1 2 3] [1 2 3] [1 2 3] [1 2 3]};

% Physical parameters
E = 2e11;
nu = 0.3;
thickness = 0.01;

problem_data.E_coeff = @(x, y, z) E * ones(size(x));
problem_data.nu_coeff = @(x, y, z) nu * ones(size(x));
problem_data.thickness = thickness;

% Source and boundary terms
hx = @(x, y, z) zeros(size(x));
hy = @(x, y, z) zeros(size(x));
hz = @(x, y, z) 80000*ones(size(x));

problem_data.f       = @(x, y, z, ind) cat(1, ...
    reshape (hx (x,y,z), [1, size(x)]), ...
    reshape (hy (x,y,z), [1, size(x)]), ...
    reshape (hz (x,y,z), [1, size(x)]));

% DISCRETIZATION PARAMETERS
clear method_data
deg = 4;
method_data.degree      = deg * [1 1];     % Degree of the splines
method_data.regularity  = (deg-2) * [1 1]; % Regularity of the splines
method_data.nsub_coarse = [4 4];           % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];           % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nquad       = (deg+1) * [1 1]; % Points for the Gaussian quadrature rule
method_data.space_type  = 'standard';
method_data.truncated   = 1;            % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag          = 'elements';
adaptivity_data.C0_est        = 1.0;
adaptivity_data.mark_param    = .25;
adaptivity_data.mark_strategy = 'GERS';
adaptivity_data.max_level     = 10;
adaptivity_data.max_ndof      = 5000;
adaptivity_data.num_max_iter  = 8;
adaptivity_data.max_nel       = 5000;
adaptivity_data.tol           = 1e-5;

% GRAPHICS
plot_data.print_info = true;
plot_data.plot_hmesh = false;

% CALL TO THE SOLVER
[geometry, hmsh, hspace, u, solution_data] = adaptivity_kirchhoff_love_shell_mp_C1 (problem_data, method_data, adaptivity_data, plot_data);

% VTK file
npts = [51 51];
output_file = 'KL_adaptivity_paraboloid.pvd';
sp_to_vtk (u, hspace, geometry, npts, output_file, 'Displacement')
