% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as multipatch AS-G1 
problem_data.geo_name = 'geo_3patch_ASG1.txt';

% Physical parameters
lambda = 5e-2;
problem_data.lambda = @(x, y) lambda* ones(size(x));
alpha = 1;
beta  = 1;
problem_data.mu  = @(x) 3 * alpha * x.^2 - beta;
problem_data.dmu = @(x) 6 * alpha * x;

% Time and time step size
problem_data.initial_time = 0;
problem_data.Time_max = 2;

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
deg = 3;
method_data.degree     = [deg deg];         % Degree of the splines
method_data.regularity = [deg-2 deg-2];     % Regularity of the splines
method_data.nquad      = [deg+1 deg+1];     % Points for the Gaussian quadrature rule
method_data.interface_regularity = 1;

% time integration parameters
method_data.rho_inf_gen_alpha = 0.5;
method_data.dt = 2e-2;

% penalty parameters
method_data.Cpen_nitsche = 1e4 * lambda; % Nitsche's method parameter
method_data.Cpen_projection = 1000;      % parameter of the penalized L2 projection (see initial conditions)

% hierarchical structure
nel = 2;
method_data.nsub_coarse = [nel nel];    % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];        % Number of subdivisions for each refinement
method_data.space_type  = 'standard';   % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 1;            % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag = 'elements';

adaptivity_data.estimator_type = 'field';
adaptivity_data.mark_param = .2;
adaptivity_data.time_delay = 0.;

adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level = 3;
adaptivity_data.max_ndof = 5000;
adaptivity_data.num_max_iter = 10;
adaptivity_data.max_nel = 500;
adaptivity_data.tol = 1e-5;

adaptivity_data.adm_type = 'T-admissible';
adaptivity_data.adm_class = 2;

% 3) INITIAL CONDITIONS 
clear initial_conditions
% initial_conditions.restart_flag = 0;
mean = 0.0;
var = 0.2;
ic_fun = @(x, y) mean + (rand(size(x))*2-1)*var;
% ic_fun = @(x, y) 0.1 * cos(2*pi*x) .* cos(2*pi*y);

problem_data.fun_u = ic_fun;
% problem_data.fun_udot = ic_fun_udot;

% create folder
folder_name = strcat('cahn_hilliard_results_p',num2str(deg),'_levels',num2str(adaptivity_data.max_level),'_lambda',num2str(lambda),'_adaptive_',adaptivity_data.estimator_type);
file_name = 'Threepatch_Cahn_Hilliard_adaptive';
status = mkdir(folder_name);

save_info.folder_name = folder_name;
save_info.file_name = file_name;
save_info.time_save = linspace(0, problem_data.Time_max, 11);

% 3) CALL TO THE SOLVER
[geometry, hmsh, hspace, results] = ...
  adaptivity_cahn_hilliard_mp_C1(problem_data, method_data, adaptivity_data, save_info);

%% 4) POST-PROCESSING

filename = strcat( folder_name,'/filenum_to_time_uniform_mesh.mat');
time_steps = results.time;
save(filename, 'time_steps');


