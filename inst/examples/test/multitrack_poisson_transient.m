% PHYSICAL DATA OF THE PROBLEM
clear all
close all
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square_x10.txt';

for iMaxLevel = 8:8
% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];

% Set non-linear flag
% if true the non-linear solver is used
problem_data.flag_nl = false;
% Set lumped matrix
problem_data.lumped = false;
% Non-linear analysis
problem_data.Newton_tol = 1.0e-06;
problem_data.num_Newton_iter = 20;
problem_data.non_linear_convergence_flag = 1;
problem_data.alpha = 1;  % alpha = 0 explicit Euler; alpha=1/2 Crank-Nicholson; aplha=1 backward Euler

% Physical parameters
problem_data.c_diff  = @(x, y) 29.0e-3 * ones(size(x));
problem_data.grad_c_diff = @(x, y) cat (1, ...
    reshape (zeros(size(x)), [1, size(x)]), ...
    reshape (zeros(size(x)), [1, size(x)]));
problem_data.c_cap = @(x, y) 650 * ones(size(x));

% Time discretization
n_tracks = 5;
n_time_steps = 500;
n_time_steps_track = n_time_steps / n_tracks;
problem_data.time_discretization = linspace(0.0, 5.0, n_time_steps);

% %Heat Source arc path
% length_x = 10;
% length_y = 10;
% theta_start = pi/4;
% theta_end = 3*pi/4;
% origin_x = 5;
% origin_y = 2.5;
% initial_radius = 5;
% hatching_space = 0.05;
% x_path = [];
% y_path = [];
% for iTrack=1:n_tracks
%     if rem(iTrack,2)==1
%         for iTimeStep=1:n_time_steps_track
%             x_path = [x_path origin_x+cos(theta_start+(iTimeStep-1)*(theta_end-theta_start)/(n_time_steps_track-1))...
%                 *(initial_radius-hatching_space*(iTrack-1))];
%             y_path = [y_path origin_y+sin(theta_start+(iTimeStep-1)*(theta_end-theta_start)/(n_time_steps_track-1))...
%                 *(initial_radius-hatching_space*(iTrack-1))];
%         end
%     else
%         for iTimeStep=1:n_time_steps_track
%             x_path = [x_path origin_x+cos(theta_end-(iTimeStep-1)*(theta_end-theta_start)/(n_time_steps_track-1))...
%                 *(initial_radius-hatching_space*(iTrack-1))];
%             y_path = [y_path origin_y+sin(theta_end-(iTimeStep-1)*(theta_end-theta_start)/(n_time_steps_track-1))...
%                 *(initial_radius-hatching_space*(iTrack-1))];
%         end
%     end
% end
% problem_data.path = [x_path', y_path'];

% Heat Source line path
length_x = 10;
length_y = 10;
origin_x = 1;
origin_y = 5-0.05*(n_tracks)/2;
laser_path_length = 8;
hatching_space = 0.050;
x_path = [];
y_path = [];
for iTrack=1:n_tracks
    if rem(iTrack,2)==1
        for iTimeStep=1:n_time_steps_track
            x_path = [x_path origin_x+(laser_path_length*(iTimeStep-1)/(n_time_steps_track-1))];
            y_path = [y_path origin_y+hatching_space*(iTrack-1)];
        end
    else
        for iTimeStep=1:n_time_steps_track
            x_path = [x_path origin_x+laser_path_length-(laser_path_length*(iTimeStep-1)/(n_time_steps_track-1))];
            y_path = [y_path origin_y+hatching_space*(iTrack-1)];
        end
    end
end
problem_data.path = [x_path', y_path'];

% Source and boundary terms
problem_data.f = gaussian_bubble_source_half_radius;
problem_data.h = @(x, y, ind) ones(size(x))*0.0;
problem_data.initial_temperature = 20;                  % [°C]

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = [3 3];                     % Degree of the splines
method_data.regularity  = method_data.degree-1;      % Regularity of the splines
method_data.nsub_coarse = [1 1];                     % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];                     % Number of subdivisions for each refinement
method_data.nquad       = method_data.degree+1;      % Points for the Gaussian quadrature rule
method_data.space_type  = 'standard';                % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 1;                         % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag = 'elements';
adaptivity_data.doCoarsening = true;
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = 0.15;
adaptivity_data.mark_param_coarsening = 0.25;
adaptivity_data.adm_strategy = 'admissible'; % 'admissible' or 'balancing'
adaptivity_data.adm = 2; % 2 if admissible 1 if balancing
adaptivity_data.radius = [1, 1];
adaptivity_data.coarsening_flag = 'any'; %'any', 'all'
adaptivity_data.coarse_flag = 'L2_global'; % 'bezier', 'MS_all', 'MS_old', 'L2_global'
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level = iMaxLevel;
adaptivity_data.max_ndof = 100000;
adaptivity_data.num_max_iter = iMaxLevel+1;
adaptivity_data.max_nel = 100000;
adaptivity_data.tol = 1.0e-01;
adaptivity_data.timeToRefine = linspace(1, n_time_steps+1, n_time_steps+1); 
   
% GRAPHICS
problem_output.folder = sprintf('geopdes_hierarchical_transient/output/Multitrack_GV15_%d', iMaxLevel);
mkdir(problem_output.folder);
plot_data.plot_hmesh = true;
plot_data.adaptivity = false;
plot_data.print_info = true;
plot_data.plot_matlab = false;
plot_data.time_steps_to_post_process = linspace(1,n_time_steps+1,n_time_steps+1);  
plot_data.file_name = strcat(problem_output.folder, '/multitrack_2D_%d');
plot_data.file_name_mesh = strcat(problem_output.folder, '/multitrack_2D_mesh_%d.png');
plot_data.file_name_results= strcat(problem_output.folder, '/multitrack_2D_results.csv');
plot_data.file_name_L2 = strcat(problem_output.folder, '/multitrack_2D_solution.csv');
plot_data.file_name_H1 = strcat(problem_output.folder, '/multitrack_2D_gradient.csv');


plot_data.npoints_x = 101;        %number of points x-direction in post-processing
plot_data.npoints_y = 101;        %number of points y-direction in post-processing
plot_data.npoints_z = 1;          %number of points z-direction in post-processing

[geometry, msh, space, u, ~] = adaptivity_poisson_transient(problem_data, method_data, adaptivity_data, plot_data);
end