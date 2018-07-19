% PHYSICAL DATA OF THE PROBLEM
clear problem_data
close all
% Physical domain, defined as NURBS map given in a text file
ndim = 2;
problem_data.geo_name = 'geo_square.txt';
mkdir('MinimalOutput');

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = 1:2*ndim;

% Physical parameters
problem_data.c_diff  = @(varargin) ones(size(varargin{1}));

% Source and boundary terms
C = 25;
problem_data.f = @(x, y) (4*C^3*(x - y)) ./ ((C*(x - y)).^2 + 1).^2;
problem_data.h = @(x, y, ind) atan (C*(x - y));

problem_data.uex = @(x, y) atan (C*(x - y));
problem_data.graduex = @(x,y) cat (1, ...
    reshape (C./(1 + (C*(x-y)).^2), [1, size(x)]), ...
    reshape (-C./(1 + (C*(x-y)).^2), [1, size(x)]));

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = 3*ones(1,ndim);       % Degree of the splines
method_data.regularity  = 2*ones(1,ndim);       % Regularity of the splines
method_data.nsub_coarse = 2*ones(1,ndim);       % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = 2*ones(1,ndim);       % Number of subdivisions for each refinement
method_data.nquad       = 4*ones(1,ndim);       % Points for the Gaussian quadrature rule
method_data.space_type  = 'standard'; % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 0;           % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag = 'elements';
% adaptivity_data.flag = 'functions';
adaptivity_data.coarsening_flag = 'all';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = .5;
adaptivity_data.coarse_flag = 'L2_global';
adaptivity_data.adm_strategy = 'admissible'; % 'admissible' or 'balancing'
adaptivity_data.adm = 1;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.mark_param_coarsening = .25;
adaptivity_data.max_level = 4;
adaptivity_data.max_ndof = 5000;
adaptivity_data.num_max_iter = 1;
adaptivity_data.max_nel = 5000;
adaptivity_data.tol = 1e-5;

% GRAPHICS
plot_data.print_info = true;
plot_data.plot_hmesh = false;
plot_data.plot_discrete_sol = false;


nlev = adaptivity_data.max_level-1;
cells = cell (nlev, 1);
cells{nlev} = 1:(2^nlev)^2;

geometry = geo_load (problem_data.geo_name);
[hmsh, hspace] = build_hspace_from_cells (ndim, method_data.degree(1), method_data.nsub_coarse(1), cells, method_data.space_type, method_data.truncated);

hm = hmsh; hs = hspace;

iter = 1;
while iter <= adaptivity_data.num_max_iter
    u = adaptivity_solve_laplace (hmsh, hspace, problem_data);
    [eh1_gal(1), el2_gal(1), eh1s_gal(1)] = sp_h1_error (hspace, hmsh, u, problem_data.uex, problem_data.graduex);
    
    % For Massimo's code
    hspace.dofs = u;
    npts = [21 21];
    [eu, F] = sp_eval (u, hspace, geometry, npts);
    %     [eu, F] = sp_eval (u, hspace, geometry, npts);
    figure(iter); surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), eu);
    title = sprintf('MinimalOutput/MinimalSolutionCoarseningFSBError20_%d.png',iter);
    saveas(gcf,title);
    fig_mesh = hmsh_plot_cells (hmsh, 20, figure(iter+100) );
    
    est = adaptivity_estimate_laplace (u, hmsh, hspace, problem_data, adaptivity_data);
    gest = norm (est);
    if (plot_data.print_info); fprintf('Computed error estimator: %f \n', gest); end
    
    [marked, num_marked] = adaptivity_mark_coarsening (est, hmsh, hspace, adaptivity_data);
    
    % COARSEN
    switch(adaptivity_data.coarse_flag)
        case 'bezier'
%             [hmsh, hspace, u] = adaptivity_coarsen (hmsh, hspace, marked, adaptivity_data);
            [hmsh, hspace, u] = adaptivity_coarsen_fsb(hmsh, hspace, marked, adaptivity_data);
        case 'MS_all'
            [hmsh, hspace, C_coarse] = adaptivity_coarsen (hmsh, hspace, marked, adaptivity_data);
            u = C_coarse * u;
        case 'MS_old'
            [hmsh, hspace, C_coarse] = adaptivity_coarsen (hmsh, hspace, marked, adaptivity_data);
            u = C_coarse * u;
        case 'L2_global'
            [hmsh, hspace, C_coarse] = adaptivity_coarsen (hmsh, hspace, marked, adaptivity_data);
            u = C_coarse * u;
    end
    
    
    [eu, F] = sp_eval (u, hspace, geometry, npts);
    %     [eu, F] = sp_eval (u, hspace, geometry, npts);
    figure(iter+20); surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), eu);
    title = sprintf('MinimalOutput/MinimalSolutionCoarseningFSBError20_%d.png',iter);
    saveas(gcf,title);
    fig_mesh = hmsh_plot_cells (hmsh, 20, figure(iter+102) );
    title = sprintf('MinimalOutput/MinimalMeshCoarseningFSBError20_%d.png',iter);
    saveas(fig_mesh,title);

    % sp_plot_solution (u_Massimo, hspace, geometry, [101 101]);
    shading interp
    
    iter = iter + 1;

end