% PHYSICAL DATA OF THE PROBLEM
close all
clear problem_data
clc

% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_3DPrintingLayer_Patil.txt';
% Generate Output folder
problem_output.folder = 'output';
mkdir(problem_output.folder);

% Set non-linear flag
% if true the non-linear solver is used
problem_data.flag_nl = true;
% Set lumped matrix
problem_data.lumped = true;

% Non-linear analysis
problem_data.Newton_tol = 1.0e-05;
problem_data.num_Newton_iter = 20;
problem_data.non_linear_convergence_flag = 1;

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];
problem_data.convection_sides = [1,2,3,4,5];
problem_data.radiation_sides = [];

% Physical parameters
problem_data.c_diff  =   @conductivity; %@(x,y,z) ones(size(x))*29;         %
problem_data.grad_c_diff =   @conductivity_der_3D; % @(x,y,z) cat (1, ...
           % reshape (zeros(size(x)), [1, size(x)]), ...
           % reshape (zeros(size(x)), [1, size(x)]), ...
           % reshape (zeros(size(x)), [1, size(x)]));
                                                            %
problem_data.c_cap = @capacity; %  @(x,y,z) ones(size(x))*7820*600;     % 
problem_data.initial_temperature = 80.0;                      %[°C]

% Time discretization
n_time_steps = 1000;
time_end = 0.0008333;
problem_data.time_discretization = linspace(0.0, time_end, n_time_steps);

% Heat Source path
laser_radius = 1.0e-04; %[m]
laser_penetration_depth = 15e-06; %[m]
absorbtion_coeff = 0.35;

problem_data.x_begin = 0.0005;
problem_data.x_end = 0.0015;
x_path_1 = linspace(problem_data.x_begin, problem_data.x_end, n_time_steps / 10);
x_path_2 = linspace(problem_data.x_end, problem_data.x_begin, n_time_steps / 10);
x_path = [x_path_1, x_path_2, x_path_1, x_path_2, x_path_1, x_path_2, x_path_1, x_path_2, x_path_1, x_path_2];

problem_data.y_begin = 0.0005;
problem_data.y_end = 0.0005;
hatch_space = laser_radius;
y_path_1 = linspace(problem_data.y_begin, problem_data.y_end, n_time_steps / 10);
y_path_2 = linspace(problem_data.y_begin + hatch_space, problem_data.y_end +hatch_space, n_time_steps / 10);
y_path_3 = linspace(problem_data.y_begin + hatch_space * 2, problem_data.y_end +hatch_space * 2, n_time_steps / 10);
y_path_4 = linspace(problem_data.y_begin + hatch_space * 3, problem_data.y_end +hatch_space * 3, n_time_steps / 10);
y_path_5 = linspace(problem_data.y_begin + hatch_space * 4, problem_data.y_end +hatch_space * 4, n_time_steps / 10);
y_path_6 = linspace(problem_data.y_begin + hatch_space * 5, problem_data.y_end +hatch_space * 5, n_time_steps / 10);
y_path_7 = linspace(problem_data.y_begin + hatch_space * 6, problem_data.y_end +hatch_space * 6, n_time_steps / 10);
y_path_8 = linspace(problem_data.y_begin + hatch_space * 7, problem_data.y_end +hatch_space * 7, n_time_steps / 10);
y_path_9 = linspace(problem_data.y_begin + hatch_space * 8, problem_data.y_end +hatch_space * 8, n_time_steps / 10);
y_path_10 = linspace(problem_data.y_begin + hatch_space * 9, problem_data.y_end +hatch_space * 9, n_time_steps / 10);

y_path = [y_path_1, y_path_2, y_path_3, y_path_4, y_path_5, y_path_6, y_path_7, y_path_8, y_path_9, y_path_10];

problem_data.z_begin = 0.000030;
problem_data.z_end = 0.000030;
z_path = linspace(problem_data.z_begin, problem_data.z_end, n_time_steps);
problem_data.path = [x_path', y_path', z_path'];



% Radiation and Convection Boundaries
problem_data.rad =  0.0;
problem_data.u_r = 3000.0;
problem_data.rad_fun = @(x,y,z,u,ind) problem_data.rad*(problem_data.u_r*ones(size(u))- u);
problem_data.rad_tilda =  0.0;
problem_data.rad_tilda_fun =  @(x,y,z,ind) ones(size(x)) * problem_data.rad_tilda;
problem_data.conv =  0.0;
problem_data.conv_fun = @(x,y,z,ind) ones(size(x)) * problem_data.conv;
problem_data.u_e = 20.0;
problem_data.conv_fun_rhs = @(x,y,z,u,ind) problem_data.conv*(problem_data.u_e*ones(size(u))- u);
problem_data.alpha = 1/2;  % alpha = 0 explicit Euler; alpha=1/2 Crank-Nicholson; aplha=1 backward Euler

% Source and boundary terms
problem_data.f = @(x,y,z,path_x,path_y,path_z) absorbtion_coeff * 6.0*sqrt(3.0)*180.0/(pi * sqrt(pi) * laser_radius * laser_radius * laser_penetration_depth) * ...
    exp(- 3*(x-path_x).^2/laser_radius^2 - 3*(y-path_y).^2/laser_radius^2 - 3*(z-path_z).^2/(laser_penetration_depth)^2);         % Body Load
problem_data.h = @(x, y, z, ind) ones(size(x))*problem_data.initial_temperature;                                                  % Dirichlet Boundaries

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = [2 2 2];        % Degree of the splines
method_data.regularity  = [1 1 1];        % Regularity of the splines
method_data.nsub_coarse = [1 1 1];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2 2];        % Number of subdivisions for each refinement
method_data.nquad       = [3 3 3];        % Points for the Gaussian quadrature rule
method_data.space_type  = 'standard';     % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 1;              % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag = 'elements';
adaptivity_data.doCoarsening = true;
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_neighbours = true;
adaptivity_data.mark_param = 0.75;
adaptivity_data.mark_param_coarsening = 0.25;
adaptivity_data.crp = 2.0;                     %coarsening relaxation parameter
adaptivity_data.radius = [laser_radius, laser_radius, laser_penetration_depth / 2]*2.0;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level = 6;
adaptivity_data.max_ndof = 100000;
adaptivity_data.num_max_iter = 5;
adaptivity_data.max_nel = 10000;
adaptivity_data.tol = 1.0e-04;
adaptivity_data.timeToRefine = linspace(1,n_time_steps,n_time_steps);

% GRAPHICS
plot_data.plot_hmesh = false;
plot_data.adaptivity = false;
plot_data.print_info = true;
plot_data.plot_matlab = true;
plot_data.time_steps_to_post_process = [1,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,...
    525,550,575,600,625,650,675,700,725,750,775,800,825,850,875,900,925,950,975,999,1000];%linspace(1,n_time_steps,n_time_steps);%[1,10,25,50,75,100,150,200,250,300,350,400];%200; % [1,25,50,75,100,150,200,250,300,350,400];% [1,5,10,15,20];
plot_data.file_name = strcat(problem_output.folder, '/Patil_FullLayer_1000_nl_maxLevel6_3D_lumped_%d.vts');
plot_data.file_name_mesh = strcat(problem_output.folder, '/Patil_FullLayer_1000_nl_maxLevel6_3D_lumped_mesh_%d');
plot_data.file_name_dofs = strcat(problem_output.folder, '/Patil_FullLayer_1000_nl_maxLevel6_3D_lumped_dofs');
plot_data.file_name_temp_plot = strcat(problem_output.folder, '/Patil_FullLayer_1000_nl_maxLevel6_3D_lumped_temperature_%d');

% plot_data.file_name_err = strcat(problem_output.folder, '/poisson_adaptivity_Fachinotti_travelling_heat_source_3D_error_%d.vts');
plot_data.npoints_x = 201;        %number of points x-direction in post-processing
plot_data.npoints_y = 101;        %number of points x-direction in post-processing
plot_data.npoints_z = 21;         %number of points x-direction in post-processing

fid = fopen (strcat(problem_output.folder, '/Patil_FullLayer_1000_nl_maxLevel6_3D_lumped_time_%d'), 'w');
tic
[geometry, hmsh, hspace, u, solution_data] = adaptivity_poisson_transient(problem_data, method_data, adaptivity_data, plot_data);
toc
fprintf (fid, num2str(toc));
fclose(fid);