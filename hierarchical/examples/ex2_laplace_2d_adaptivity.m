% PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));

% Source and boundary terms
C = 100;
normax2 = @(x,y) ((x-.5).^2+(y-.5).^2);
problem_data.uex = @(x,y) exp(-C*normax2(x,y));
problem_data.f = @(x,y) 4*C*(1-C*normax2(x,y)).*problem_data.uex(x,y);
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) problem_data.uex(x,y);

% Exact solution (optional)
problem_data.uex =@(x,y) exp(-C*normax2(x,y));
problem_data.graduex = @(x,y) -2*C*cat (1, ...
            reshape (problem_data.uex(x,y).*(x-.5), [1, size(x)]), ...
            reshape (problem_data.uex(x,y).*(y-.5), [1, size(x)]));
        

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = [3 3];       % Degree of the splines
method_data.regularity  = [2 2];       % Regularity of the splines
method_data.nsub_coarse = [2 2];       % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];       % Number of subdivisions for each refinement
method_data.nquad       = [4 4];       % Points for the Gaussian quadrature rule
method_data.space_type  = 0;           % 0: , 1: Full basis (B-splines)
method_data.truncated   = 0;           % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag = 'elements';
%adaptivity_data.flag = 'functions';
adaptivity_data.mark_param = .5;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level = 10;
adaptivity_data.max_ndof = 5000;
adaptivity_data.num_max_iter = 15;
adaptivity_data.max_nel = 5000;
adaptivity_data.tol = 1e-5;

% GRAPHICS
plot_hmesh = true;
plot_discrete_sol = true;

[geometry, hmsh, hspace, u, gest, err_h1s, iter] = adaptivity_solve_laplace(problem_data, method_data, adaptivity_data, plot_hmesh, plot_discrete_sol);
