clear problem_data method_data

nrb(1) = nrbdegelev(nrb4surf([-1 -1], [0 -1], [-1 0], [0 0]), [1 1]);
nrb(2) = nrbdegelev(nrb4surf([0 -1], [1 -1], [0 0], [1 0]), [1 1]);
nrb(3) = nrbdegelev(nrb4surf([-1 0], [0 0], [-1 1], [0 1]), [1 1]);
nrb(4) = nrbdegelev(nrb4surf([0 0], [1 0], [0 1], [1 1]), [1 1]);
nrb(1).coefs(3,:,:) = [-1 0 0; 0 1 1; 0 1 1];
nrb(2).coefs(3,:,:) = [0 1 1; 0 1 1; -1 0 0];
nrb(3).coefs(3,:,:) = [0 0 -1; 1 1 0; 1 1 0];
nrb(4).coefs(3,:,:) = [1 1 0; 1 1 0; 0 0 -1];

problem_data.geo_name = nrb;

problem_data.drchlt_sides = [1 2 3 4 5 6 7 8];
problem_data.drchlt_components = {[1 2 3] [1 2 3] [1 2 3] [1 2 3] [1 2 3] [1 2 3] [1 2 3] [1 2 3]};

% Physical parameters
E = 2e11;
nu = 0.3;
thickness = 0.01;

problem_data.E_coeff = @(x, y, z) E * ones(size(x));
problem_data.nu_coeff = @(x, y, z) nu * ones(size(x));
problem_data.thickness = thickness;

% Source and boundary terms
hx = @(x, y, z) zeros(size(x));
hy = @(x, y, z) zeros(size(x));
hz = @(x, y, z) 80000*ones(size(x));

problem_data.f       = @(x, y, z, ind) cat(1, ...
    reshape (hx (x,y,z), [1, size(x)]), ...
    reshape (hy (x,y,z), [1, size(x)]), ...
    reshape (hz (x,y,z), [1, size(x)]));

% Discretization parameters
method_data.degree      = [5 5];        % Degree of the splines
method_data.regularity  = [3 3];        % Regularity of the splines
method_data.nsub_coarse = [4 4];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nquad       = method_data.degree+1;       % Points for the Gaussian quadrature rule
method_data.space_type  = 'standard';
method_data.truncated   = 1;            % 0: False, 1: True

adaptivity_data.flag          = 'elements';
adaptivity_data.C0_est        = 1.0;
adaptivity_data.mark_param    = .5;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level     = 10;
adaptivity_data.max_ndof      = 5000;
adaptivity_data.num_max_iter  = 8;
adaptivity_data.max_nel       = 5000;
adaptivity_data.tol           = 1e-5;

% GRAPHICS
plot_data.print_info = true;
plot_data.plot_hmesh = false;
plot_data.plot_discrete_sol = false;

% CALL TO THE SOLVER
[geometry, hmsh, hspace, u, solution_data] = adaptivity_kirchhoff_love_shell_mp_C1 (problem_data, method_data, adaptivity_data, plot_data);

% % VTK file
% output_file = sprintf('output_vtk/paraboloid_%ipatch', numel(nrb));
% vtk_pts = {linspace(0, 1, 51), linspace(0, 1, 51)};
% fprintf ('The result is saved in the file %s \n \n', output_file);
% 
% lambda_lame = @(x, y, z) ((problem_data.nu_coeff(x,y,z).*problem_data.E_coeff(x,y,z)) ./ ...
%   ((1+problem_data.nu_coeff(x,y,z)).*(1-2*problem_data.nu_coeff(x,y,z)))); 
% mu_lame = @(x, y, z) (problem_data.E_coeff(x,y,z)/(2*(1+problem_data.nu_coeff(x,y,z))) * ones (size (x)));
% sp_to_vtk (u, space, geometry, vtk_pts, output_file, {'Displacement', 'Gradient', 'Stress'}, {'value', 'gradient', 'stress'}, lambda_lame, mu_lame)

% % Evaluate solution in one point
% [Cpatch, cols] = sp_compute_Cpatch_vector (space, 1, msh.rdim);
% u_loc = Cpatch * u(cols);
% sp_vec = sp_vector (repmat (space.sp_patch(1), 1, msh.rdim), msh.msh_patch{1});
% eu2{ii} = sp_eval(u_loc,sp_vec,geometry(1), {1 1});
