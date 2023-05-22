clc
close all
clear all





% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
% nrb1 = nrb4surf ([0 0], [.5 0], [0 .5], [.5 .5]);
% nrb2 = nrb4surf ([.5 0], [1 0], [.5 .5], [1 .5]);
% nrb3 = nrb4surf ([0 0.5], [.5 0.5], [0 1], [.5 1]);
% nrb4 = nrb4surf ([.5 0.5], [1 0.5], [.5 1], [1 1]);
% nrb(1) = nrb1;
% nrb(2) = nrb2;
% nrb(3) = nrb3;
% nrb(4) = nrb4;

problem_data.geo_name = 'curved_3patch_bicubic.txt' ;% nrb; %'geo_square_mp.txt';

% Physical parameters
lambda = 5e-2;%(1/(4*sqrt(2) *pi))^2; %6.15e-4;% 2.45e-3; 
problem_data.lambda = @(x, y) lambda* ones(size(x));


% Time and time step size
Time_max = 50;
dt = 1e-1;
problem_data.time = 0;
problem_data.Time_max = Time_max;


% penalty parameters
problem_data.Cpen_nitsche = 1e4 * lambda; % Nitsche's method parameter
problem_data.Cpen_projection = 1000;      % parameter of the penalized L2 projection (see initial conditions)


% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
nel =1;
p= 3;
method_data.degree     = [p p];         % Degree of the splines
method_data.regularity = [p-2 p-2];     % Regularity of the splines
%method_data.nsub       = [nel nel];     % Number of subdivisions
method_data.nquad      = [p+1 p+1];     % Points for the Gaussian quadrature rule

% time integration parameters
method_data.rho_inf_gen_alpha = 0.5;
method_data.dt =dt;

% hierarchical structure
method_data.nsub_coarse = [nel nel];    % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];        % Number of subdivisions for each refinement
method_data.space_type  = 'standard';   % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 1;            % 0: False, 1: True


% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag = 'elements';
adaptivity_data.C0_est = 1.0;


adaptivity_data.estimator_type = 'field';
adaptivity_data.mark_param = .1;
adaptivity_data.mark_param_coarsening = .1;
adaptivity_data.adm_class = p;
adaptivity_data.time_delay = 0.;


adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level = 5;
adaptivity_data.max_ndof = 5000;
adaptivity_data.num_max_iter = 10;
adaptivity_data.max_nel = 500;
adaptivity_data.tol = 1e-5;


adaptivity_data.adm_type = 'T-admissible';

% 3) INITIAL CONDITIONS 
clear initial_conditions
mean = 0.4;
var = 0.005;
ic_fun = @(x, y) mean + (rand(size(x))*2-1)*var;
%ic_fun = @(x, y) 0.1 * cos(2*pi*x) .* cos(2*pi*y);


%ic_fun = load("initial_conditions_nel_64_64_p_2_2.mat");
%ic_fun = ic_fun.u_0;
initial_conditions.fun_u = ic_fun;


% create folder
folder_name = strcat('chan_hilliard_results_p',num2str(p),'_nel',num2str(nel),'_lambda',num2str(lambda),'_adaptive_',adaptivity_data.estimator_type);
status = mkdir(folder_name);
status = rmdir(folder_name);
status = mkdir(folder_name);

save_info.folder_name = folder_name;
save_info.time_save = linspace(-.000001,problem_data.Time_max,51);

% 3) CALL TO THE SOLVER
[geometry, hmsh, hspace, results] = adaptivity_cahn_hilliard_mp_C1(problem_data, method_data, adaptivity_data, initial_conditions, save_info);


%% 4) POST-PROCESSING



filename = strcat( folder_name,'/filenum_to_time_uniform_mesh.mat');
time_steps = results.time;
save(filename, 'time_steps');




    
    













%% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW
% u = results.u;
% output_file = 'cahn_hilliard_adaptive';
% 



% % 5.2) EXPORT TO PARAVIEW    
% vtk_pts = {linspace(0, 1, 40), linspace(0, 1, 40)};
% for step = 1:length(results.time)
%   output_file = strcat( folder_name,'/Square_cahn_hilliard_', num2str(step) );
%   fprintf ('The result is saved in the file %s \n \n', output_file);
%   sp_to_vtk (results.u(:,step), hspace, geometry, vtk_pts, output_file, {'u','grad_u'}, {'value','gradient'})
% end
%     

% npts = [51 51];
% fprintf ('The result is saved in the file %s \n \n', output_file);
% sp_to_vtk (u, hspace, geometry, npts, output_file,{'u', 'grad_u'}, {'value', 'gradient'})
% 
% 
% [eu, F] = sp_eval (u, hspace, geometry, npts, {'value', 'gradient'});
% fig=figure; subplot (1,3,1)
% surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), eu{1}, 'EdgeColor', 'none')
% view(0,90);
% colorbar
% title('solution')
% shading interp
% subplot(1,3,2)
% surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), squeeze ( eu{2}(1,:,:) ), 'EdgeColor', 'none');
% view(0,90);
% title('grad x')
% shading interp
% subplot(1,3,3)
% surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), squeeze ( eu{2}(2,:,:) ), 'EdgeColor', 'none');
% view(0,90);
% title('grad y')
% shading interp
% saveas(fig , 'solution_hierarchical.png')
% %%
% fig=figure;
% value = sqrt(eu{2}(1,:,:).^2 + eu{2}(2,:,:).^2);
% %value = value/max(value,[],"all");
% surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), squeeze ( value ), 'EdgeColor', 'none', 'FaceAlpha',0.5);
% hold on
% hmsh_plot_cells (hmsh)
% view(0,90);
% colorbar
% shading interp
% axis equal tight
% saveas(fig , 'indicator.png')
% 
% 
% fig=figure;
% hmsh_plot_cells (hmsh)
% view(0,90);
% shading interp
% axis equal tight
% saveas(fig , 'mesh_hierarchical.png')



%!demo
%! ex_laplace_square

%!test
%! problem_data.geo_name = 'geo_square.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4];
%! problem_data.c_diff  = @(x, y) ones(size(x));
%! problem_data.f = @(x, y) zeros (size (x));
%! problem_data.g = @test_square_g_nmnn;
%! problem_data.h = @(x, y, ind) exp (x) .* sin(y);
%! problem_data.uex     = @(x, y) exp (x) .* sin (y);
%! problem_data.graduex = @(x, y) cat (1, ...
%!                       reshape (exp(x).*sin(y), [1, size(x)]), ...
%!                       reshape (exp(x).*cos(y), [1, size(x)]));
%! method_data.degree     = [3 3];       % Degree of the splines
%! method_data.regularity = [2 2];       % Regularity of the splines
%! method_data.nsub       = [9 9];       % Number of subdivisions
%! method_data.nquad      = [4 4];       % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u] = solve_laplace (problem_data, method_data);
%! [error_h1, error_l2] = ...
%!           sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex);
%! assert (msh.nel, 81)
%! assert (space.ndof, 144)
%! assert (error_h1, 9.86428525677199e-06, 1e-14)
%! assert (error_l2, 1.68004134750130e-07, 1e-14)