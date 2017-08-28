% EX_PLANE_STRAIN_ADAPTIVITY: solve the plane-strain problem on one quarter of a cylinder.

% PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = nrb4surf([0 0], [1 0], [0 1], [1 1]);

% Type of boundary conditions
% problem_data.nmnn_sides   = [2];
% problem_data.drchlt_sides = [];
% problem_data.press_sides  = [1];
% problem_data.symm_sides   = [3 4];
problem_data.nmnn_sides   = [];
problem_data.press_sides  = [];
problem_data.drchlt_sides = [1 2 3 4];
problem_data.symm_sides   = [];

% Physical parameters
E  =  1; nu = .3; 
%E  =  1; nu = 0; 
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
% P = 1;
% problem_data.f = @(x, y) zeros (2, size (x, 1), size (x, 2));
% problem_data.g = @(x, y, ind) test_plane_strain_ring_g_nmnn (x, y, P, nu, ind);
% problem_data.h = @(x, y, ind) test_plane_strain_ring_uex (x, y, E, nu, P);
% problem_data.p = @(x, y, ind) P * ones (size (x));
fx = @(x, y) -(-(problem_data.mu_lame(x, y)*3 + problem_data.lambda_lame(x, y)).*sin(2*pi*x).*sin(2*pi*y) + ...
     (problem_data.mu_lame(x, y) + problem_data.lambda_lame(x, y)).*cos(2*pi*x).*cos(2*pi*y))*(2*pi)^2;
fy = fx;
problem_data.f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));
problem_data.h       = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));

% Exact solution (optional)
%problem_data.uex = @(x, y) test_plane_strain_ring_uex (x, y, E, nu, P);
uxex = @(x,y) sin(2*pi*x).*(sin(2*pi*y));
uyex = @(x,y) sin(2*pi*x).*(sin(2*pi*y));
problem_data.uex = @(x, y) cat(1, ...
                reshape (uxex (x,y), [1, size(x)]), ...
                reshape (uyex (x,y), [1, size(x)]));
problem_data.h = @(x,y,ind) problem_data.uex(x,y);

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree      = [3 3];        % Degree of the splines
method_data.regularity  = [2 2];        % Regularity of the splines
method_data.nsub_coarse = [2 2];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];        % Number of subdivisions for each refinement
method_data.nquad       = [4 4];        % Points for the Gaussian quadrature rule
method_data.space_type  = 'simplified'; % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 0;            % 0: False, 1: True

adaptivity_data.flag = 'elements';
% adaptivity_data.flag = 'functions';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = .5;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level = 10;
adaptivity_data.max_ndof = 15000;
adaptivity_data.num_max_iter = 8;
adaptivity_data.max_nel = 5000;
adaptivity_data.tol = 1e-5;

% GRAPHICS
plot_data.print_info = true;
plot_data.plot_hmesh = false;
plot_data.plot_discrete_sol = false;

% 3) CALL TO THE SOLVER
[geometry, hmsh, hspace, u] = adaptivity_linear_elasticity (problem_data, method_data, adaptivity_data, plot_data);

% 4) POST-PROCESSING. 
% 4.1) Export to Paraview
output_file = 'plane_strain_ring_Deg3_Reg2_Sub9';

vtk_pts = {linspace(0, 1, 21), linspace(0, 1, 21)};
fprintf ('results being saved in: %s \n \n', output_file)
sp_to_vtk (u, hspace, geometry, vtk_pts, output_file, {'displacement', 'stress'}, {'value', 'stress'}, ...
    problem_data.lambda_lame, problem_data.mu_lame)

% 4.2) Plot in Matlab. Comparison with the exact solution.
[eu, F] = sp_eval (u, hspace, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

subplot (1,2,1)
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
title ('Numerical solution'), axis equal tight
subplot (1,2,2)
eu2 = problem_data.uex (X, Y);
quiver (X, Y, squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
title ('Exact solution'), axis equal tight

error_l2 = sp_l2_error (hspace, hmsh, u, problem_data.uex)

