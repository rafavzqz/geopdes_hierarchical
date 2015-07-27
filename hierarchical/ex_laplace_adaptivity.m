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
problem_data.f = @(x, y) zeros (size (x));
problem_data.g = @test_square_g_nmnn;
problem_data.h = @(x, y, ind) exp (x) .* sin(y);

% Exact solution (optional)
problem_data.uex     = @(x, y) exp (x) .* sin (y);
problem_data.graduex = @(x, y) cat (1, ...
                       reshape (exp(x).*sin(y), [1, size(x)]), ...
                       reshape (exp(x).*cos(y), [1, size(x)]));

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree     = [3 3];       % Degree of the splines
method_data.regularity = [2 2];       % Regularity of the splines
method_data.nsub       = [2 2];       % Number of subdivisions
method_data.nquad      = [4 4];       % Points for the Gaussian quadrature rule
method_data.space_type = 0;           % 0: , 1: Full basis (B-splines)
method_data.truncated = 0;            % 0: False, 1: True

% ADAPTIVITY PARAMETERS

adaptivity_data.est_type = 3;
adaptivity_data.mark_param = .5;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level = 10;
adaptivity_data.max_ndof = 5000;
adaptivity_data.num_max_iter = 15;
adaptivity_data.max_nel = 5000;
