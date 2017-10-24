% PHYSICAL DATA OF THE PROBLEM
close all
clear problem_data
clc
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_2DPrintingLayer.txt';
% Generate Output folder
problem_output.folder = 'output';
mkdir(problem_output.folder);

% Set non-linear flag
% if true the non-linear solver is used
problem_data.flag_nl = true;
problem_data.lumped = false;
problem_data.non_linear_convergence_flag = 1;

% Non-linear analysis
problem_data.Newton_tol = 1.0e-05;
problem_data.num_Newton_iter = 20;
problem_data.alpha = 1/2;  % alpha = 0 explicit Euler; alpha=1/2 Crank-Nicholson; aplha=1 backward Euler

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];
problem_data.convection_sides = [1,2,3,4];
problem_data.radiation_sides = 1;

% Physical parameters
problem_data.c_diff  =   @conductivity; %    
problem_data.grad_c_diff = @conductivity_der;%  @(x,y) cat (1,reshape (zeros(size(x)), [1, size(x)]), reshape (zeros(size(x)), [1, size(x)])); % 
problem_data.c_cap =  @capacity;% 
problem_data.initial_temperature = 20; %[°C]


% Time discretization
n_time_steps = 30;
time_end = 20.0;
problem_data.time_discretization = linspace(0.0, time_end, n_time_steps + 1);

% Radiation and Convection Boundaries
problem_data.rad =  0.15;
problem_data.u_r = 3000.0;
problem_data.rad_fun = @(x,y,u,ind) problem_data.rad*(problem_data.u_r*ones(size(u))- u);
problem_data.rad_tilda =  0.0;
problem_data.rad_tilda_fun =  @(x,y,ind) ones(size(x)) * problem_data.rad_tilda;
problem_data.conv =  2.0;
problem_data.conv_fun = @(x,y,ind) ones(size(x)) * problem_data.conv;
problem_data.u_e = 20.0;
problem_data.conv_fun_rhs = @(x,y,u,ind) problem_data.conv*(problem_data.u_e*ones(size(u))- u);

% Heat Source path
x_begin = 0.05;
x_end = 0.15;
x_path = linspace(x_begin, x_end, n_time_steps+1);
y_begin = 0.1;
y_end = 0.1;
y_path = linspace(y_begin, y_end, n_time_steps+1);
problem_data.path = [x_path', y_path'];

% Source and boundary terms
problem_data.f = moving_heat_source;        % Body Load
problem_data.h = [];                        % Dirichlet Boundaries
problem_data.g = [];                        % Neumann Boundaries

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = [3 3];        % Degree of the splines
method_data.regularity  = [2 2];        % Regularity of the splines
method_data.nsub_coarse = [2 2];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];        % Number of subdivisions for each refinement
method_data.nquad       = method_data.degree+1;        % Points for the Gaussian quadrature rule
method_data.space_type  = 'standard';                  % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 1;                           % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag = 'elements';
adaptivity_data.C0_est = 1.0;
adaptivity_data.doCoarsening = false;
adaptivity_data.mark_param = 0.75;
adaptivity_data.mark_param_coarsening = 0.25;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.radius = [0.01, 0.01]*2.0;
adaptivity_data.adm = 2;
adaptivity_data.max_level = 5;
adaptivity_data.max_ndof = 15000;
adaptivity_data.num_max_iter = 5;
adaptivity_data.max_nel = 15000;
adaptivity_data.tol = 5.0e-01;
adaptivity_data.timeToRefine = linspace(1, n_time_steps, n_time_steps); 

% GRAPHICS
plot_data.plot_hmesh = true;
plot_data.adaptivity = false;
plot_data.plot_discrete_sol = false;
plot_data.print_info = true;
plot_data.plot_matlab = false;
plot_data.time_steps_to_post_process =linspace(1,n_time_steps,n_time_steps); %  [1, 5, 15, 25, 35, 45, 50, 100, 150, 200, 250, 300, 350, 400]; %        [10, 20]; %
plot_data.file_name = strcat(problem_output.folder, '/travelling_heat_source_nonlinear_2D_%d.vts');
plot_data.file_name_mesh = strcat(problem_output.folder, '/travelling_heat_source_nonlinear_2D_%d');
plot_data.file_name_dofs = strcat(problem_output.folder, '/travelling_heat_source_nonlinear_2D_dofs');
plot_data.npoints_x = 101;        %number of points x-direction in post-processing
plot_data.npoints_y = 101;        %number of points y-direction in post-processing
plot_data.npoints_z = 1;          %number of points z-direction in post-processing

[geometry, hmsh, hspace, u, solution_data] = adaptivity_poisson_transient(problem_data, method_data, adaptivity_data, plot_data);