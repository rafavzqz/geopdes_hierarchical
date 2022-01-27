% PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file

% 8 patches (L-shape)
p1=[0 0]; p2=[1 0]; p3=[2/3 1/3]; p4=[-1/3 1/3]; p5=[0 -1]; p6=[-1/3 -2/3]; p7=[-1 -1]; p8=[-2/3 -2/3]; p9=[-1 1]; p10=[-2/3 2/3]; p11=[1 1]; p12=[2/3 2/3];
nrb(1) = nrb4surf(p1,p2,p4,p3);
nrb(2) = nrb4surf(p1,p5,p4,p6);
nrb(3) = nrb4surf(p6,p5,p8,p7);
nrb(4) = nrb4surf(p7, p9, p8,p10);
nrb(5) = nrb4surf(p8,p10,p6,p4);
nrb(6) = nrb4surf(p10,p12,p4,p3);
nrb(7) = nrb4surf(p9,p11,p10,p12);
nrb(8) = nrb4surf(p3,p12,p2,p11);
problem_data.geo_name = nrb;

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4 5 6];
problem_data.weak_drchlt_sides = [];
        
% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));
% problem_data.grad_c_diff = @(x, y) cat (1, ...  
%             reshape (zeros(size(x)), [1, size(x)]), ...
%             reshape (zeros(size(x)), [1, size(x)]));

% Source and boundary terms
C = 1;
p = -1;
% Bilaplacian
problem_data.f = @bilaplacian_rhs_Lshape2; %zeros(size(x));%@bilaplacian_rhs_Lshape2; %@(x, y) p*ones(size(x)); %@(x, y) zeros(size(x)); %bilaplacian
problem_data.g = @(x, y, ind) bilaplacian_Lshape_g_nmnn_r_8patches(x,y,ind); %@(x, y,iside) zeros(size(x)); % @bilaplacian_Lshape_g_nmnn_r2; %@bilaplacian_Lshape_g_nmnn_r2_C; %@bilaplacian_Lshape_g_nmnn_r;%pass here the normal derivative (on the boundary)
problem_data.h = @(x, y, ind) solution_bilaplacian_Lshape (x, y);%(x, y, ind) ones(size(x)); %@(x, y, ind) solution_bilaplacian_Lshape (x, y); %@(x, y,iside) zeros(size(x));  %solution on the boundary

% Exact solution (optional)
problem_data.uex     = @(x, y) solution_bilaplacian_Lshape (x, y);
problem_data.graduex = @(x, y) solution_bilaplacian_Lshape_grad (x, y);
problem_data.hessuex = @(x, y) solution_bilaplacian_Lshape_hessian (x, y);


% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = [3 3];        % Degree of the splines
method_data.regularity  = [1 1];        % Regularity of the splines
method_data.nsub_coarse = [4 4];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];        % Number of subdivisions for each refinement
method_data.nquad      = method_data.degree+1;       % Points for the Gaussian quadrature rule
method_data.Cpen = 10 * (min(method_data.degree) + 1);
method_data.space_type  = 'standard'; % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 1;            % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag = 'elements';
adaptivity_data.C0_est = 10;
adaptivity_data.mark_param = .25;
adaptivity_data.mark_strategy = 'GERS';
adaptivity_data.max_level = 6;
adaptivity_data.max_ndof = 15000;
adaptivity_data.num_max_iter = 8;
adaptivity_data.max_nel = 15000;
adaptivity_data.tol = 1e-7;
adaptivity_data.adm_class = 2;
adaptivity_data.adm_type = 'T-admissible';

% GRAPHICS
plot_data.print_info = true;
plot_data.plot_hmesh = false;
plot_data.plot_discrete_sol = false;

[geometry, hmsh, hspace, uu, solution_data] = adaptivity_bilaplace_mp_C1 (problem_data, method_data, adaptivity_data, plot_data);

% save solution_data_bilaplacian solution_data adaptivity_data method_data problem_data
% hspace.C_L2 = [];
% save solution_bilaplacian hmsh hspace

%% % EXPORT VTK FILE
%subplot(1,2,1)
%npts = [51 51];
%sp_plot_solution (uu, hspace, geometry, hmsh, npts); shading interp
%vtk_pts = {linspace(0,1,npts(1)), linspace(0,1,npts(2))};
%subplot(1,2,2)
%for iptc = 1:hmsh.npatch
%  F = reshape (geometry(iptc).map(vtk_pts), [2, npts]);
%  X = reshape (F(1,:), npts); Y = reshape (F(2,:), npts);
%  surf(X, Y, problem_data.uex(X,Y));
%  hold on; shading interp
%end
% 
% warning ('on','geopdes:nrbmultipatch')
% 
% %MESH PLOT
% figure(2)
% hmsh_plot_cells (hmsh, 2, 2);
% figure (2)
% view(0,90)
% set(gcf,'color','white')
% axis equal
% %axis off
