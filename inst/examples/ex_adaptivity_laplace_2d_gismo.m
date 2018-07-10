% PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_rectangle.txt';
filename = '/Users/ondine/Documents/gismo/filedata/domain2d/squareTHB.xml';
problem_data_gismo.geo_name = gsTHBSpline(filename);

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];
problem_data_gismo.nmnn_sides = problem_data.nmnn_sides;
problem_data_gismo.drchlt_sides = problem_data.drchlt_sides;

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));
problem_data.grad_c_diff = @(x, y) cat (1, ...
            reshape (zeros(size(x)), [1, size(x)]), ...
            reshape (zeros(size(x)), [1, size(x)]));
problem_data_gismo.c_diff = problem_data.c_diff;
problem_data_gismo.grad_c_diff = problem_data.grad_c_diff;

% Source and boundary terms
problem_data.f = @(x,y) 2*(2*pi)^2*sin(2*pi*x).*sin(2*pi*y);
problem_data.h = @(x,y, ind) zeros (size (x));
problem_data_gismo.f = problem_data.f;
problem_data_gismo.h = problem_data.h;

% Exact solution (optional)
problem_data.uex     = @(x,y) sin(2*pi*x).*sin(2*pi*y);
problem_data.graduex = @(x, y) cat (1, ...
                       reshape (2*pi*cos(2*pi*x).*sin(2*pi*y), [1, size(x)]), ...
                       reshape (2*pi*sin(2*pi*x).*cos(2*pi*y), [1, size(x)]));
problem_data_gismo.uex = problem_data.uex;
problem_data_gismo.graduex = problem_data.graduex;

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = [2, 2];            % Degree of the splines
method_data.regularity  = [1, 1];            % Regularity of the splines
method_data.nsub_coarse = [4, 4];            % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2, 2];            % Number of subdivisions for each refinement
method_data.nquad       = [3, 3];            % Points for the Gaussian quadrature rule
method_data.space_type  = 'standard';        % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 1;                 % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag = 'elements';
% adaptivity_data.flag = 'functions';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = .5;
adaptivity_data.mark_strategy = 'MS'; %GERS
adaptivity_data.max_level = 12;
adaptivity_data.max_ndof = 100000;
adaptivity_data.num_max_iter = 15;
adaptivity_data.max_nel = 100000; 
adaptivity_data.tol = 1e-3;

% GRAPHICS
plot_data.plot_hmesh = false;
plot_data.plot_discrete_sol = false;
plot_data.print_info = true;

[geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);
[geometry_gismo, hmsh_gismo, hspace_gismo, u_gismo, solution_data_gismo] = adaptivity_laplace (problem_data_gismo, method_data, adaptivity_data, plot_data);

% EXPORT VTK FILE
npts = [51 51];
output_file = 'laplace_adaptivity_square_geopdes.vts';
output_file_gismo = 'laplace_adaptivity_square_gismo.vts';
sp_to_vtk (u, hspace, geometry, npts, output_file, {'solution', 'gradient'}, {'value', 'gradient'})
sp_to_vtk (u_gismo, hspace_gismo, geometry_gismo, npts, output_file_gismo, {'solution', 'gradient'}, {'value', 'gradient'})

% Plot in Octave/Matlab
[eu, F] = sp_eval (u, hspace, geometry, npts);
figure; subplot (1,2,1)
surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), eu)
subplot(1,2,2)
surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), squeeze (problem_data.uex(F(1,:,:), F(2,:,:))));

[eu_gismo, F_gismo] = sp_eval (u_gismo, hspace_gismo, geometry_gismo, npts);
figure; subplot (1,2,1)
surf (squeeze(F_gismo(1,:,:)), squeeze(F_gismo(2,:,:)), eu_gismo)
subplot(1,2,2)
surf (squeeze(F_gismo(1,:,:)), squeeze(F_gismo(2,:,:)), squeeze (problem_data_gismo.uex(F_gismo(1,:,:), F_gismo(2,:,:))));

assert(length(u) == length(u_gismo))
[eu_err, F_err] = sp_eval (u-u_gismo, hspace, geometry, npts);
figure; subplot (1,2,1)
surf (squeeze(F_err(1,:,:)), squeeze(F_err(2,:,:)), eu_err)
subplot(1,2,2)
surf (squeeze(F_err(1,:,:)), squeeze(F_err(2,:,:)), squeeze (problem_data.uex(F_err(1,:,:), F_err(2,:,:))));
