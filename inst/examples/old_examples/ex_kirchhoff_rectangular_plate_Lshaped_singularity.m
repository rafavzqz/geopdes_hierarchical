
% PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Geometry definition (unrefined geometry)
base = 1.0; height = 1.0;
p11 =[0 0]; p12 =[base 0]; p21 =[-1.0 height]; p22 =[base height];
srf = nrb4surf (p11,p12,p21,p22);

problem_data.geo_name = srf;

% Boundary conditions
problem_data.simply_supported_sides = [1 3];
problem_data.clamped_sides = [];
problem_data.nmnn_sides = [];
problem_data.press_sides = [];
problem_data.drchlt_sides = [2 4];
problem_data.prescribed_moment_sides = [2 4];

problem_data.D = @(x, y) ones(size(x));

% Singular function
problem_data.f = @(x, y) zeros (size (x));
problem_data.h = @(x, y, ind) singular_function_bilaplace (x, y);
problem_data.M  = @(x, y, iside) analyticalPrescribedMoment (x, y, iside);
% Exact solution (optional)
problem_data.uex     = @(x, y) singular_function_bilaplace (x, y);
problem_data.graduex = @(x, y) singular_function_bilaplace_gradient (x, y);
problem_data.hessuex = @(x, y) singular_function_bilaplace_hessian (x, y);

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
adaptivity_data.mark_param = .8;
adaptivity_data.mark_strategy = 'MS';
% adaptivity_data.mark_strategy = 'GR';
adaptivity_data.flag = 'elements';
% adaptivity_data.flag = 'functions';
adaptivity_data.num_max_iter = 17; %21
adaptivity_data.max_nel = 5000;
adaptivity_data.tol = 1e-10;        

problem_data.geometry = geo_load (problem_data.geo_name);
[knots, zeta] = kntrefine (problem_data.geometry.nurbs.knots, method_data.nsub_coarse-1, method_data.degree, method_data.regularity);

rule     = msh_gauss_nodes (method_data.nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh   = msh_cartesian (zeta, qn, qw, problem_data.geometry, 'der2', true);
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
% output_file = sprintf('~/Local/Src/GeoPDs/Kirchhoff_Plate_Simply_Supported_Singular_Adaptivity_%d',iter );
% vtk_pts = {linspace(0, 1, 31), linspace(0, 1, 31)};
% fprintf ('The result is saved in the file %s \n \n', output_file);
% sp_to_vtk (u, hspace, problem_data.geometry, vtk_pts, output_file, 'u')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% ESTIMATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  estimator = adaptivity_estimate_bilaplace (u, hmsh, hspace, problem_data, adaptivity_data);
  if (isfield (problem_data, 'hessuex'))
%     [err_h2(iter), ~] = sp_l2_error (hspace, hmsh, u, problem_data.uex);
    [~, ~, ~, err_h2(iter), ~, ~, ~, ~, ~, ~] = sp_h2_error (hspace, hmsh, u, ...
                                            problem_data.uex, problem_data.graduex, problem_data.hessuex);
    dof_vector(iter) = ndof;
  end
%   if (isfield (problem_data, 'graduex'))
%     [err_h1(iter), ~, ~] = sp_h1_error (hspace, hmsh, u, problem_data.uex, problem_data.graduex);
%     dof_vector(iter) = ndof;
%   end

%   space_bubble = space_bubble_function (hmsh, 'plate');
%   estimator = compute_bubble_estimator_bilaplacian (hspace, space_bubble, hmsh, problem_data.f, [], [], u, problem_data.D);
  
  gest(iter) = norm (estimator);
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% MARK
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% output_file = sprintf('~/Local/Src/GeoPDs/Kirchhoff_Plate_Corner_Singularity_Adaptivity_%d',iter+1 );
% vtk_pts = {linspace(0, 1, 31), linspace(0, 1, 31)};
% fprintf ('The result is saved in the file %s \n \n', output_file);
% sp_to_vtk (u, hspace, problem_data.geometry, vtk_pts, output_file, 'u')
% hmsh_plot_cells(hmsh);

estimator = adaptivity_estimate_bilaplace (u, hmsh, hspace, problem_data, adaptivity_data);

% space_bubble = space_bubble_function (hmsh, 'plate');
% estimator = compute_bubble_estimator_bilaplacian (hspace, space_bubble, hmsh, problem_data.f, [], [], u, problem_data.D);

gest(iter+1) = norm (estimator);

if (isfield (problem_data, 'hessuex'))
    [~, ~, ~, err_h2(iter+1), ~, ~, ~, ~, ~, ~] = sp_h2_error (hspace, hmsh, u, ...
                                        problem_data.uex, problem_data.graduex, problem_data.hessuex);
    dof_vector(iter+1) = hspace.ndof;
end
% if (isfield (problem_data, 'graduex'))
%     [err_h1(iter), ~, ~] = sp_h1_error (hspace, hmsh, u, problem_data.uex, problem_data.graduex);
%     dof_vector(iter+1) = hspace.ndof;
% end

data = [dof_vector.^0.5; err_h2; gest];
fileID = fopen('data_corner_singularity_2x2_bubble.txt','w');
fprintf(fileID,'%12.14f\t %12.14f\t %12.14f\n',data);
fclose(fileID);

disp('END:')

function res = analyticalPrescribedMoment (x, y, iside)
    [theta, rad] = cart2pol (x, y);
    theta = (theta < 0).*(2*acos(-1) + theta) + (theta >= 0) .* theta;
    switch iside
      case 1
        disp('Not implemented');
      case 2
        res = -8/9*rad.^(-2/3).*cos(theta).*sin(theta).*cos(4/3*theta) + ...
            4/9*rad.^(-2/3).*sin(4/3*theta).*(cos(theta).^2 - sin(theta).^2);
      case 3
        res = 8/9*rad.^(-2/3).*cos(theta).*sin(theta).*cos(4/3*theta) - ...
            4/9*rad.^(-2/3).*sin(4/3*theta).*(cos(theta).^2 - sin(theta).^2);
      case 4
        res = 8/9*rad.^(-2/3).*cos(theta).*sin(theta).*cos(4/3*theta) - ...
            4/9*rad.^(-2/3).*sin(4/3*theta).*(cos(theta).^2 - sin(theta).^2);  
      otherwise
        disp('Unknown case')
    end
end


function f = singular_function_bilaplace (x, y)
  [theta, rad] = cart2pol (x, y);
  theta = (theta < 0).*(2*acos(-1) + theta) + (theta >= 0) .* theta;
  f = rad.^(4/3) .* sin((4/3)*theta);
end

function f = singular_function_bilaplace_gradient (x, y)
  f = zeros(2, size(x,1), size(x,2));
  [theta, rad] = cart2pol (x, y);
  theta = (theta < 0).*(2*acos(-1) + theta) + (theta >= 0) .* theta;
  f(1,:,:) = (4/3)*rad.^(1/3) .* sin((1/3)*theta);
  f(2,:,:) = (4/3)*rad.^(1/3) .* cos((1/3)*theta);
end

function f = singular_function_bilaplace_hessian (x, y)
  f = zeros(4, size(x,1), size(x,2));
  [theta, rad] = cart2pol (x, y);
  theta = (theta < 0).*(2*acos(-1) + theta) + (theta >= 0) .* theta;
  f(1,:,:) = -8/9*rad.^(-2/3).*cos(theta).*sin(theta).*cos(4/3*theta) + ...
            4/9*rad.^(-2/3).*sin(4/3*theta).*(cos(theta).^2 - sin(theta).^2);
  f(2,:,:) = 8/9*rad.^(-2/3).*cos(theta).*sin(theta).*sin(4/3*theta) + ...
            4/9*rad.^(-2/3).*cos(4/3*theta).*(cos(theta).^2 - sin(theta).^2);
  f(3,:,:) = 8/9*rad.^(-2/3).*cos(theta).*sin(theta).*sin(4/3*theta) + ...
            4/9*rad.^(-2/3).*cos(4/3*theta).*(cos(theta).^2 - sin(theta).^2);
  f(4,:,:) = 8/9*rad.^(-2/3).*cos(theta).*sin(theta).*cos(4/3*theta) - ...
            4/9*rad.^(-2/3).*sin(4/3*theta).*(cos(theta).^2 - sin(theta).^2);
end