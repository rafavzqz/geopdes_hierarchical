clear problem_data
close all
clc
% Initial domain
problem_data.geo_name = nrbline ([0 0], [1.0 0]);

% CHOICE OF THE DISCRETIZATION PARAMETERS (Initial mesh)
clear method_data
method_data.degree      = 2;                % Degree of the splines
method_data.regularity  = 1;                % Regularity of the splines
method_data.nsub_coarse = 3;                % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = 2;                % Number of subdivisions for each refinement
method_data.nquad       = 3;                % Points for the Gaussian quadrature rule
method_data.space_type  = 'standard';       % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 1;                % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag = 'elements';
% adaptivity_data.flag = 'functions';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = 0.75;
adaptivity_data.mark_param_coarsening = 0.01;
adaptivity_data.adm_strategy = 'balancing'; % 'admissible' or 'balancing'
adaptivity_data.adm = 1;
adaptivity_data.coarse_flag = 'bezier';
adaptivity_data.max_level = 5;
adaptivity_data.max_ndof = 15000;
adaptivity_data.num_max_iter = 20;
adaptivity_data.max_nel = 15000;
adaptivity_data.tol = 1.0e-03;

% GRAPHICS
plot_data.plot_hmesh = false;
plot_data.adaptivity = false;
plot_data.plot_discrete_sol = false;
plot_data.print_info = false;
plot_data.plot_matlab = false;
plot_data.npoints_x = 100;       %number of points x-direction in post-processing
plot_data.npoints_y = 1;         %number of points y-direction in post-processing
plot_data.npoints_z = 1;         %number of points z-direction in post-processing

%% Generate the first level of initial mesh
[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (problem_data.geo_name);

[knots, zeta] = kntrefine (geometry.nurbs.knots, method_data.nsub_coarse-1, method_data.degree, method_data.regularity);
% tensor product level 1
rule     = msh_gauss_nodes (method_data.nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh   = msh_cartesian (zeta, qn, qw, geometry);
space = sp_bspline (knots, method_data.degree, msh);
% hierarchical space
hmsh     = hierarchical_mesh (msh, method_data.nsub_refine);
hspace   = hierarchical_space (hmsh, space, method_data.space_type, method_data.truncated);

hspace.dofs = [1 2 3 4 5]';
initial_values = hspace.dofs;


%% Add second and third level using THB-refinement
% Second level
marked_ref = cell(1, hspace.nlevels);
marked_ref{1} = [1 2];
marked_ref{2} = [];
[hmsh, hspace, C_ref] = adaptivity_refine (hmsh, hspace, marked_ref, adaptivity_data);

% hmsh_plot_cells (hmsh, 20, (figure(1)) );

hspace.dofs = C_ref*hspace.dofs;
temp_dofs = hspace.dofs;

% add new level

hmsh.nlevels = hmsh.nlevels + 1;
hmsh.active{hmsh.nlevels} = [];
hmsh.deactivated{hmsh.nlevels} = [];
hmsh.nel_per_level(hmsh.nlevels) = 0;
hmsh.mesh_of_level(hmsh.nlevels) = msh_refine (hmsh.mesh_of_level(hmsh.nlevels-1), hmsh.nsub);
hmsh.msh_lev{hmsh.nlevels} = [];

% Third level
marked_ref = cell(1, hspace.nlevels);
marked_ref{1} = [];
marked_ref{2} = [1 2 3];
marked_ref{3} = [];
[hmsh, hspace, C_ref] = adaptivity_refine (hmsh, hspace, marked_ref, adaptivity_data);

% hmsh_plot_cells (hmsh, 20, (figure(2)) );

hspace.dofs = C_ref*hspace.dofs;
temp_dofs2 = hspace.dofs;


% add new level
hmsh.nlevels = hmsh.nlevels + 1;
hmsh.active{hmsh.nlevels} = [];
hmsh.deactivated{hmsh.nlevels} = [];
hmsh.nel_per_level(hmsh.nlevels) = 0;
hmsh.mesh_of_level(hmsh.nlevels) = msh_refine (hmsh.mesh_of_level(hmsh.nlevels-1), hmsh.nsub);
hmsh.msh_lev{hmsh.nlevels} = [];


% % plot initial state
% npts = [plot_data.npoints_x];
% [eu, F] = sp_eval (hspace.dofs, hspace, geometry, npts);
% figure(4); plot (squeeze(F(1,:,:)), eu)

%% Refinement =============================================================
marked_ref{1} = [];
marked_ref{2} = [];
marked_ref{3} = [3 4 5 6];
marked_ref{4} = [];
[hmsh, hspace, Cref] = adaptivity_refine (hmsh, hspace, marked_ref, adaptivity_data);

% hmsh_plot_cells (hmsh, 20, (figure(3)));
hspace.dofs = Cref*hspace.dofs;

% % plot refined state
% npts = [plot_data.npoints_x];
% [eu, F] = sp_eval (hspace.dofs, hspace, geometry, npts);
% figure(5); plot (squeeze(F(1,:,:)), eu)


%% Coarsening back to the initial state ===================================
marked_coarse{1} = [];
marked_coarse{2} = [];
marked_coarse{3} = [];
marked_coarse{4} = [5 6 7 8];
[hmsh, hspace, u] = adaptivity_coarsen(hmsh, hspace, marked_coarse, adaptivity_data);

% hmsh_plot_cells (hmsh, 20, (figure(3)));
hspace.dofs = u;

knots = hspace.space_of_level(hmsh.nlevels).knots;
length = size(hspace.Csub{hmsh.nlevels},1);
cps = [linspace(1,30,length); linspace(1,30,length)*0.2; zeros(1,length);...
    ones(1,length)];
curv = nrbmak(cps,knots{1});

R   = RefinementOperator(crv,hspace);

ML  = MultiLevelOperator(hspace,R);

MLC = MultiLevelExtractionOperator(crv,hspace,ML);

neval_points = 100;
nrbplot(curv, neval_points);
eval_points = linspace(0,1,neval_points);
[B, id] = nrbbasisfun (eval_points, curv);
bsplines = zeros(size(hspace.Csub{hmsh.nlevels},1),neval_points);
for ip = 1:neval_points
    bsplines(id(ip,:),ip)=transpose(B(ip,:));
end
plot(eval_points, bsplines);

% [ C ] = bzrextr( knots{1}, curv.order-1 );
% bernstein = zeros((numel(knots{1})-2*(curv.order-1))*curv.order,neval_points);
% for ip = 1:neval_points
%     bernstein(1+(id(ip,1)-1)*curv.order:(id(ip,1)-1)*curv.order+curv.order,ip)= C(:,:,id(ip,1))\B(ip,:)';
% end
% plot(eval_points, bernstein);

THB = MLC{hmsh.nlevels}*bsplines;
plot(eval_points, THB);

marked_coarse{1} = [];
marked_coarse{2} = [];
marked_coarse{3} = [1 2 3 4];
marked_coarse{4} = [];
[hmsh, hspace, u] = adaptivity_coarsen(hmsh, hspace, marked_coarse, adaptivity_data);
hmsh_plot_cells (hmsh, 20, (figure(4)));
hspace.dofs = u;

knots = hspace.space_of_level(hmsh.nlevels).knots;
cps = [(hspace.Csub{hmsh.nlevels}*hspace.dofs)'; (hspace.Csub{hmsh.nlevels}*hspace.dofs)'*0.2; zeros(1, size(hspace.Csub{hmsh.nlevels},1));...
    ones(1, size(hspace.Csub{hmsh.nlevels},1))];
curv = nrbmak(cps,knots{1});
neval_points = 1000;
nrbplot(curv, neval_points);
eval_points = linspace(0,1,neval_points);
[B, id] = nrbbasisfun (eval_points, curv);
bsplines = zeros(size(hspace.Csub{hmsh.nlevels},1),neval_points);
for ip = 1:neval_points
    bsplines(id(ip,:),ip)=transpose(B(ip,:));
end
plot(eval_points, bsplines);
THB = hspace.Csub{hmsh.nlevels}'*bsplines;
plot(eval_points, THB);

marked_coarse{1} = [];
marked_coarse{2} = [1 2];
marked_coarse{3} = [];
marked_coarse{4} = [];
[hmsh, hspace, u] = adaptivity_coarsen(hmsh, hspace, marked_coarse, adaptivity_data);
hmsh_plot_cells (hmsh, 20, (figure(5)));
hspace.dofs = u;

knots = hspace.space_of_level(hmsh.nlevels).knots;
cps = [(hspace.Csub{hmsh.nlevels}*hspace.dofs)'; (hspace.Csub{hmsh.nlevels}*hspace.dofs)'*0.2; zeros(1, size(hspace.Csub{hmsh.nlevels},1));...
    ones(1, size(hspace.Csub{hmsh.nlevels},1))];
curv = nrbmak(cps,knots{1});
neval_points = 20;
nrbplot(curv, neval_points);
eval_points = linspace(0,1,neval_points);
[B, id] = nrbbasisfun (eval_points, curv);
bsplines = zeros(size(hspace.Csub{hmsh.nlevels},1),neval_points);
for ip = 1:neval_points
    bsplines(id(ip,:),ip)=transpose(B(ip,:));
end
plot(eval_points, bsplines);
THB = hspace.Csub{hmsh.nlevels}'*bsplines;
plot(eval_points, THB);

% plot coarse state
npts = [plot_data.npoints_x];
[eu, F] = sp_eval (u, hspace, geometry, npts);
figure(6); plot (squeeze(F(1,:,:)), eu)
