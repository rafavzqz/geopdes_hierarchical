% EX_KIRCHHOFF_RECTANGULAR_PLATE_CLAMPED_ADAPTIVITY: solve bilaplacian in a clamped rectangular plate.
% PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Geometry definition (unrefined geometry)
base = 1.0; height = 1.0;
p11 =[0 0]; p12 =[base 0]; p21 =[0 height]; p22 =[base height];
srf = nrb4surf (p11,p12,p21,p22);

problem_data.geo_name = srf;

a = 5;
b = 5;

c = 0; %14
d = 0;
% Boundary conditions
problem_data.simply_supported_sides = [];
problem_data.clamped_sides = [];
problem_data.nmnn_sides = [];
problem_data.press_sides = [];
% problem_data.drchlt_sides = [];
% problem_data.prescribed_moment_sides = [];
problem_data.drchlt_sides = [1 2 3 4];
problem_data.prescribed_moment_sides = [1 2 3 4];

problem_data.h = @(x, y, ind)  x.^a.*y.^b.*(x - 1).^c.*(y - 1).^d;
problem_data.M  = @(x, y, iside) analyticalPrescribedMoment (x, y, a, b, c, d, iside);

% Physical parameters
problem_data.D = @(x, y) ones(size(x));
% Source term  

problem_data.f = @(x, y) a*x.^(a - 4).*y.^b.*(a - 1).*(a - 2).*(a - 3).*(x - 1).^c.*(y - 1).^d + 4*a*c*x.^(a - 3).*y.^b.*(a - 1).*(a - 2).*(x - 1).^(c - 1).*(y - 1).^d + ...
                         6*a*c*x.^(a - 2).*y.^b.*(a - 1).*(c - 1).*(x - 1).^(c - 2).*(y - 1).^d + 4*a*c*x.^(a - 1).*y.^b.*(c - 1).*(c - 2).*(x - 1).^(c - 3).*(y - 1).^d + ...
                         c*x.^a.*y.^b.*(c - 1).*(c - 2).*(c - 3).*(x - 1).^(c - 4).*(y - 1).^d + ...
                         2*(a*b*x.^(a - 2).*y.^(b - 2).*(a - 1).*(b - 1).*(x - 1).^c.*(y - 1).^d + 2*a*b*c*x.^(a - 1).*y.^(b - 2).*(b - 1).*(x - 1).^(c - 1).*(y - 1).^d + ... 
                         2*a*b*d*x.^(a - 2).*y.^(b - 1).*(a - 1).*(x - 1).^c.*(y - 1).^(d - 1) + 4*a*b*c*d*x.^(a - 1).*y.^(b - 1).*(x - 1).^(c - 1).*(y - 1).^(d - 1) + ...
                         a*d*x.^(a - 2).*y.^b.*(a - 1).*(d - 1).*(x - 1).^c.*(y - 1).^(d - 2) + b*c*x.^a.*y.^(b - 2).*(b - 1).*(c - 1).*(x - 1).^(c - 2).*(y - 1).^d + ...
                         2*a*c*d*x.^(a - 1).*y.^b.*(d - 1).*(x - 1).^(c - 1).*(y - 1).^(d - 2) + 2*b*c*d*x.^a.*y.^(b - 1).*(c - 1).*(x - 1).^(c - 2).*(y - 1).^(d - 1) + ...
                         c*d*x.^a.*y.^b.*(c - 1).*(d - 1).*(x - 1).^(c - 2).*(y - 1).^(d - 2) ) + ...
                         b*x.^a.*y.^(b - 4).*(b - 1).*(b - 2).*(b - 3).*(x - 1).^c.*(y - 1).^d + 4*b*d*x.^a.*y.^(b - 3).*(b - 1).*(b - 2).*(x - 1).^c.*(y - 1).^(d - 1) + ...
                         6*b*d*x.^a.*y.^(b - 2).*(b - 1).*(d - 1).*(x - 1).^c.*(y - 1).^(d - 2) + 4*b*d*x.^a.*y.^(b - 1).*(d - 1).*(d - 2).*(x - 1).^c.*(y - 1).^(d - 3) + ...
                         d*x.^a.*y.^b.*(d - 1).*(d - 2).*(d - 3).*(x - 1).^c.*(y - 1).^(d - 4);

problem_data.uex     = @(x, y) x.^a.*y.^b.*(x - 1).^c.*(y - 1).^d;

problem_data.graduex = @(x, y) cat (1, ...
                       reshape (a*x.^(a - 1).*y.^b.*(x - 1).^c.*(y - 1).^d + c*x.^a.*y.^b.*(x - 1).^(c - 1).*(y - 1).^d, [1, size(x)]), ...
                       reshape (b*x.^a.*y.^(b - 1).*(x - 1).^c.*(y - 1).^d + d*x.^a.*y.^b.*(x - 1).^c.*(y - 1).^(d - 1), [1, size(x)]));

problem_data.hessuex = @(x, y) cat (1, ...
                       reshape (a*x.^(a - 2).*y.^b.*(a - 1).*(x - 1).^c.*(y - 1).^d + 2*a*c*x.^(a - 1).*y.^b.*(x - 1).^(c - 1).*(y - 1).^d + ...
                       c*x.^a.*y.^b.*(c - 1).*(x - 1).^(c - 2).*(y - 1).^d, [1, size(x)]), ...
                       reshape (a*b*x.^(a - 1).*y.^(b - 1).*(x - 1).^c.*(y - 1).^d + a*d*x.^(a - 1).*y.^b.*(x - 1).^c.*(y - 1).^(d - 1) + ...
                       b*c*x.^a.*y.^(b - 1).*(x - 1).^(c - 1).*(y - 1).^d + c*d*x.^a.*y.^b.*(x - 1).^(c - 1).*(y - 1).^(d - 1), [1, size(x)]), ...
                       reshape (a*b*x.^(a - 1).*y.^(b - 1).*(x - 1).^c.*(y - 1).^d + a*d*x.^(a - 1).*y.^b.*(x - 1).^c.*(y - 1).^(d - 1) + ...
                       b*c*x.^a.*y.^(b - 1).*(x - 1).^(c - 1).*(y - 1).^d + c*d*x.^a.*y.^b.*(x - 1).^(c - 1).*(y - 1).^(d - 1), [1, size(x)]), ...
                       reshape (b*x.^a.*y.^(b - 2).*(b - 1).*(x - 1).^c.*(y - 1).^d + 2*b*d*x.^a.*y.^(b - 1).*(x - 1).^c.*(y - 1).^(d - 1) + ...
                       d*x.^a.*y.^b.*(d - 1).*(x - 1).^c.*(y - 1).^(d - 2), [1, size(x)]));

% problem_data.uex     = @(x, y) x.^a.*y.^b;
% 
% problem_data.graduex = @(x, y) cat (1, ...
%                        reshape (a*x.^(a - 1).*y.^b, [1, size(x)]), ...
%                        reshape (b*x.^a.*y.^(b - 1), [1, size(x)]));
% 
% problem_data.hessuex = @(x, y) cat (1, ...
%                        reshape (a*x.^(a - 2).*y.^b.*(a - 1), [1, size(x)]), ...
%                        reshape (a*b*x.^(a - 1).*y.^(b - 1), [1, size(x)]), ...
%                        reshape (a*b*x.^(a - 1).*y.^(b - 1), [1, size(x)]), ...
%                        reshape (b*x.^a.*y.^(b - 2).*(b - 1), [1, size(x)]));
% 
% problem_data.f = @(x, y) a*x.^(a - 4).*y.^b.*(a - 1).*(a - 2).*(a - 3) + ...
%                          b*x.^a.*y.^(b - 4).*(b - 1).*(b - 2).*(b - 3) + ...
%                          2*(a*b*x.^(a - 2).*y.^(b - 2).*(a - 1).*(b - 1));


n = 2;
deg = 4;

clear method_data
method_data.degree      = [deg deg];        % Degree of the splines
method_data.regularity  = [deg-1 deg-1];        % Regularity of the splines
method_data.nsub_coarse = [n n];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];        % Number of subdivisions for each refinement
method_data.nquad       = [deg+1 deg+1];        % Points for the Gaussian quadrature rule
method_data.space_type  = 'simplified'; % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 0;            % 0: False, 1: True

clear adaptivity_data
adaptivity_data.mark_param = .6;
adaptivity_data.mark_strategy = 'MS';
% adaptivity_data.mark_strategy = 'GR';
adaptivity_data.flag = 'elements';
% adaptivity_data.flag = 'functions';
adaptivity_data.num_max_iter = 5; %RB 22, B 44
adaptivity_data.max_nel = 5000;
adaptivity_data.tol = 1e-10;        

problem_data.geometry = geo_load (problem_data.geo_name);
[knots, zeta] = kntrefine (problem_data.geometry.nurbs.knots, method_data.nsub_coarse-1, method_data.degree, method_data.regularity);

rule     = msh_gauss_nodes (method_data.nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh   = msh_cartesian (zeta, qn, qw, problem_data.geometry, 'der2', true);
space = sp_bspline (knots, method_data.degree, msh);
hmsh     = hierarchical_mesh (msh, method_data.nsub_refine);
hspace   = hierarchical_space (hmsh, space, method_data.space_type, method_data.truncated);

nel = zeros (1, adaptivity_data.num_max_iter+1); ndof = nel; gest = nel+NaN;
dof_vector = zeros(1, adaptivity_data.num_max_iter+1);
if (isfield (problem_data, 'hessuex'))
    err_h2 = gest;
end

%% ADAPTIVE LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter = 0;
while (iter < adaptivity_data.num_max_iter)
  iter = iter + 1;
  
  ndof = hspace.ndof;
  
  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',iter);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% SOLVE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('SOLVE:')
  % Assemblying stiffness matrix
  % CALL TO THE SOLVER
  u = adaptivity_solve_bilaplace_gradgrad_2d_iso(hspace, hmsh, problem_data);
  
% POST-PROCESSING
% EXPORT TO PARAVIEW
% output_file = sprintf('~/Local/Src/GeoPDs/Kirchhoff_Plate_Simply_Supported_Singular_Adaptivity_%d',iter );
% vtk_pts = {linspace(0, 1, 31), linspace(0, 1, 31)};
% fprintf ('The result is saved in the file %s \n \n', output_file);
% sp_to_vtk (u, hspace, problem_data.geometry, vtk_pts, output_file, 'u')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% ESTIMATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  estimator = adaptivity_estimate_bilaplace (u, hmsh, hspace, problem_data, adaptivity_data);
  if (isfield (problem_data, 'hessuex'))
%     [err_h2(iter), ~] = sp_l2_error (hspace, hmsh, u, problem_data.uex);
    [~, ~, ~, err_h2(iter), ~, ~, ~, ~, ~, ~] = sp_h2_error (hspace, hmsh, u, ...
                                            problem_data.uex, problem_data.graduex, problem_data.hessuex);
    dof_vector(iter) = ndof;
  end

%   space_bubble = space_bubble_function (hmsh, 'plate');
%   estimator = compute_bubble_estimator_bilaplacian (hspace, space_bubble, hmsh, problem_data.f, [], [], u, problem_data.D);
  
  gest(iter) = norm (estimator);
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% MARK
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  disp('MARK:')
    [marked, num_marked] = adaptivity_mark (estimator, hmsh, hspace, adaptivity_data);
  fprintf('%d %s marked for refinement \n', num_marked, adaptivity_data.flag);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% REFINE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('REFINE:')
    [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
  fprintf('\n');
  
end

u = adaptivity_solve_bilaplace_gradgrad_2d_iso(hspace, hmsh, problem_data);
% output_file = sprintf('~/Local/Src/GeoPDs/Kirchhoff_Plate_Simply_Supported_Singular_Adaptivity_%d',iter+1 );
% vtk_pts = {linspace(0, 1, 31), linspace(0, 1, 31)};
% fprintf ('The result is saved in the file %s \n \n', output_file);
% sp_to_vtk (u, hspace, problem_data.geometry, vtk_pts, output_file, 'u')
% hmsh_plot_cells(hmsh);

estimator = adaptivity_estimate_bilaplace (u, hmsh, hspace, problem_data, adaptivity_data);

% space_bubble = space_bubble_function (hmsh, 'plate');
% estimator = compute_bubble_estimator_bilaplacian (hspace, space_bubble, hmsh, problem_data.f, [], [], u, problem_data.D);

gest(iter+1) = norm (estimator);

if (isfield (problem_data, 'hessuex'))
%     [err_h2(iter), ~] = sp_l2_error (hspace, hmsh, u, problem_data.uex);
    [~, ~, ~, err_h2(iter+1), ~, ~, ~, ~, ~, ~] = sp_h2_error (hspace, hmsh, u, ...
                                        problem_data.uex, problem_data.graduex, problem_data.hessuex);
    dof_vector(iter+1) = hspace.ndof;
end

data = [dof_vector.^0.5; err_h2; gest];
fileID = fopen('data.txt','w');
fprintf(fileID,'%12.14f\t %12.14f\t %12.14f\n',data);
fclose(fileID);

disp('END:')

function res = analyticalPrescribedMoment (x, y, a, b, c, d, iside)
%     switch iside
%       case 1
%         res = a*x.^(a - 2).*y.^b.*(a - 1);
%       case 2
%         res = a*x.^(a - 2).*y.^b.*(a - 1);
%       case 3
%         res = b*x.^a.*y.^(b - 2).*(b - 1);
%       case 4
%         res = b*x.^a.*y.^(b - 2).*(b - 1);  
%       otherwise
%         disp('Unknown case')
%     end
%     switch iside
%       case 1
%         res = a*x.^(a - 2).*y.^b.*(a - 1).*(x - 1).^c.*(y - 1).^d + 2*a*c*x.^(a - 1).*y.^b.*(x - 1).^(c - 1).*(y - 1).^d + ...
%                        c*x.^a.*y.^b.*(c - 1).*(x - 1).^(c - 2).*(y - 1).^d;
%       case 2
%         res = a*x.^(a - 2).*y.^b.*(a - 1).*(x - 1).^c.*(y - 1).^d + 2*a*c*x.^(a - 1).*y.^b.*(x - 1).^(c - 1).*(y - 1).^d + ...
%                        c*x.^a.*y.^b.*(c - 1).*(x - 1).^(c - 2).*(y - 1).^d;
%       case 3
%         res = b*x.^a.*y.^(b - 2).*(b - 1).*(x - 1).^c.*(y - 1).^d + 2*b*d*x.^a.*y.^(b - 1).*(x - 1).^c.*(y - 1).^(d - 1) + ...
%                        d*x.^a.*y.^b.*(d - 1).*(x - 1).^c.*(y - 1).^(d - 2);
%       case 4
%         res = b*x.^a.*y.^(b - 2).*(b - 1).*(x - 1).^c.*(y - 1).^d + 2*b*d*x.^a.*y.^(b - 1).*(x - 1).^c.*(y - 1).^(d - 1) + ...
%                        d*x.^a.*y.^b.*(d - 1).*(x - 1).^c.*(y - 1).^(d - 2);  
%       otherwise
%         disp('Unknown case')
%     end
%     
     switch iside
      case 1
        res = 20*x.^3.*y.^5;
      case 2
        res = 20*x.^3.*y.^5;
      case 3
        res = 20*y.^3.*x.^5;
      case 4
        res = 20*y.^3.*x.^5;  
      otherwise
        disp('Unknown case')
    end
   
end