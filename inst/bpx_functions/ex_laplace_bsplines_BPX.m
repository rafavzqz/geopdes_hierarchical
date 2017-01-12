% EX_LAPLACE_SQUARE: solve the Poisson problem in the unit square with a B-spline discretization.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file

for dim = 2;

if (dim == 1)
  problem_data.geo_name = nrbline ([0 0], [1 0]);
elseif (dim == 2)
  problem_data.geo_name = nrb4surf ([0 0], [1 0], [0 1], [1 1]);
elseif (dim == 3)
  problem_data.geo_name = nrbextrude (nrb4surf ([0 0], [1 0], [0 1], [1 1]), [0 0 1]);
end

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = 1:(2*dim);

% Physical parameters
% % % % Source and boundary terms
% % % problem_data.f = @(x, y) zeros (size (x));
% % % problem_data.g = @test_square_g_nmnn;
% % % problem_data.h = @(x, y, ind) exp (x) .* sin(y);


% Not working (non-homogeneous BC)
if (dim == 2)
  k = 1; % Constant that characterizes the singularity
  problem_data.f = @(x, y) zeros(size(x));
  problem_data.h = @(x, y, ind) singular_function_laplace (x, y, k);

% Exact solution (optional)
  problem_data.uex     = @(x, y) singular_function_laplace (x, y, k);
  problem_data.graduex = @(x, y) cat (1, ...
                       reshape (-(2*y.*cos((2*atan(y./x))/3) - 2*x.*sin((2*atan(y./x))/3))./(3*(x.^2 + y.^2).^(2/3)), [1, size(x)]), ...
                       reshape ((2*x.*cos((2*atan(y./x))/3) + 2*y.*sin((2*atan(y./x))/3))./(3*(x.^2 + y.^2).^(2/3)), [1, size(x)]));
else
  problem_data.f = @(varargin) ones (size(varargin{1}));
  problem_data.h = @(varargin) zeros (size(varargin{1}));
  problem_data.uex     = @(varargin) zeros (size(varargin{1}));
  problem_data.graduex = @(varargin) zeros ([dim, size(varargin{1})]);
end
 
            
% Not working (non-homogeneous BC)
%problem_data.uex     = @(x, y) singular_function_laplace (x, y, k);
%problem_data.graduex = @(x, y) singular_function (x, y, k);

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = 3 * ones (dim, 1);       % Degree of the splines
method_data.regularity = 2 * ones (dim, 1);       % Regularity of the splines
method_data.nsub       = 1 * ones (dim, 1);       % Number of subdivisions of coarse geometry
method_data.nquad      = 4 * ones (dim, 1);       % Points for the Gaussian quadrature rule

% method_data.mesh_type = 'dyadic_non_constant';
% method_data.nint_coarse = 4; % Number of intervals of the zero level mesh
method_data.nlevels = 10; % It will also compute all the intermediate cases

method_data.mesh_type = 'uniform';

for deg = 2:4
   method_data.degree = deg * ones(dim, 1);
   method_data.regularity = method_data.degree - 1;
   
   method_data.nquad = method_data.degree + 1;
   problem_data.filename = ['results_bsplines_dim',num2str(dim),'_degree', num2str(deg)];
% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_laplace_bsplines_BPX_new (problem_data, method_data);

end

end
% % 4) POST-PROCESSING
% % 4.1) EXPORT TO PARAVIEW
% 
% output_file = 'Square_tsplines';
% 
% vtk_pts = {linspace(0, 1, 40), linspace(0, 1, 40)};
% fprintf ('The result is saved in the file %s \n \n', output_file);
% sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')
% 
% % 4.2) PLOT IN MATLAB. COMPARISON WITH THE EXACT SOLUTION
% 
% [eu, F] = sp_eval (u, space, geometry, vtk_pts);
% [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
% subplot (1,2,1)
% surf (X, Y, eu)
% title ('Numerical solution'), axis tight
% subplot (1,2,2)
% surf (X, Y, problem_data.uex (X,Y))
% title ('Exact solution'), axis tight
% 
% % Display errors of the computed solution in the L2 and H1 norm
% [error_h1, error_l2] = ...
%            sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)
% 
% %!demo
% %! ex_laplace_square
