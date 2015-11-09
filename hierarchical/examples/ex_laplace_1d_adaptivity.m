% PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = nrbline ([0 0], [1 0]);

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2];

% Physical parameters
problem_data.c_diff  = @(x) ones(size(x));

% Source and boundary terms
problem_data.f = @(x) (2*pi)^2*sin(2*pi*x);
problem_data.h = @(x, ind) zeros (size (x));

% Exact solution (optional)
problem_data.uex     = @(x) sin(2*pi*x);
problem_data.graduex = @(x) 2*pi*cos(2*pi*x);

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = 2;            % Degree of the splines
method_data.regularity  = 1;            % Regularity of the splines
method_data.nsub_coarse = 2;            % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = 2;            % Number of subdivisions for each refinement
method_data.nquad       = 3;            % Points for the Gaussian quadrature rule
method_data.space_type  = 'simplified'; % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 0;            % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
%adaptivity_data.flag = 'elements';
adaptivity_data.flag = 'functions';
adaptivity_data.mark_param = .015;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level = 12;
adaptivity_data.max_ndof = 5000;
adaptivity_data.num_max_iter = 15;
adaptivity_data.max_nel = 5000;
adaptivity_data.tol = 1e-9;

% GRAPHICS
plot_hmesh = false;
plot_discrete_sol = false;

[geometry, hmsh, hspace, u, gest, err_h1s, iter] = ...
    adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_hmesh, plot_discrete_sol);
