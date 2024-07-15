nel =1;

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square.txt';

% Physical parameters
lambda = 6.15e-4;
problem_data.lambda = @(x, y) lambda* ones(size(x));

% time  
problem_data.initial_time = 0;
problem_data.Time_max = 0.021;    

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
p = 2;
method_data.degree     = [p p];         % Degree of the splines
method_data.regularity = [p-1 p-1];     % Regularity of the splines
%method_data.nsub       = [nel nel];     % Number of subdivisions
method_data.nquad      = [p+1 p+1];     % Points for the Gaussian quadrature rule

% time integration parameters
method_data.rho_inf_gen_alpha = 0.5;
method_data.dt = 1e-3;

% penalty parameters
method_data.pen_nitsche = 1e4 * lambda; % Nitsche's  penalty constant 
method_data.pen_projection = 1000;      % penalty value to constrain the flux in L2 projection

% hierarchical structure
nel = 1;
method_data.nsub_coarse = [nel nel];    % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];        % Number of subdivisions for each refinement
method_data.space_type  = 'standard';   % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 1;            % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag = 'elements';
adaptivity_data.C0_est = 1.0;

adaptivity_data.estimator_type = 'field';
adaptivity_data.mark_param = .2;
adaptivity_data.mark_param_coarsening = .2;
adaptivity_data.time_delay = 0.;

adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level = 5;
adaptivity_data.max_ndof = 5000;
adaptivity_data.num_max_iter = 10;
adaptivity_data.max_nel = 500;
adaptivity_data.tol = 1e-5;
adaptivity_data.adm_type = 'T-admissible';
adaptivity_data.adm_class = 2;

% 3) INITIAL CONDITIONS 
clear initial_conditions
mean = 0.4;
var = 0.005;
ic_fun = @(x, y) mean + (rand(size(x))*2-1)*var;
ic_fun = @(x, y) 0.1 * cos(2*pi*x) .* cos(2*pi*y);


%ic_fun = load("initial_conditions_nel_64_64_p_2_2.mat");
%ic_fun = ic_fun.u_0;
initial_conditions.fun_u = ic_fun;


% create folder
folder_name = strcat('chan_hilliard_results_p',num2str(p),'_nel',num2str(nel),'_lambda',num2str(lambda),'_adaptive_',adaptivity_data.estimator_type);
status = mkdir(folder_name);
status = rmdir(folder_name);
status = mkdir(folder_name);

save_info.folder_name = folder_name;
save_info.time_save = linspace(-.000001,problem_data.Time_max,3);

% 3) CALL TO THE SOLVER
[geometry, hmsh, hspace, results] = adaptivity_cahn_hilliard(problem_data, method_data, adaptivity_data, initial_conditions, save_info);


%% 4) POST-PROCESSING
filename = strcat( folder_name,'/filenum_to_time_uniform_mesh.mat');
time_steps = results.time;
save(filename, 'time_steps');

