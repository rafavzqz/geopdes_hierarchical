% Adaptive isogeometric method for d-dimensional problems (d = 1,2,3)

clear hmsh hspace u est
% close all
% add_my_paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PHYSICAL DATA OF THE PROBLEM
clear problem_data  
problem = 6; 
problem_data = load_problem_data(problem);

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = [3 3];       % Degree of the splines
method_data.regularity  = method_data.degree-1;       % Regularity of the splines
method_data.nsub_coarse = [2 2];       % Number of subdivisions
method_data.nsub_refine = [2 2];       % Number of subdivisions
method_data.nquad       = method_data.degree+1;       % Points for the Gaussian quadrature rule
method_data.space_type  = 'simplified'; % 'simplified' (only children functions) or 'classical' (full basis)

% ADAPTIVITY PARAMETERS

adaptivity_data.est_type = 3;
adaptivity_data.mark_param = .5;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level = 10;
adaptivity_data.max_ndof = 5000;
adaptivity_data.num_max_iter = 15;
adaptivity_data.max_nel = 5000;
adaptivity_data.tol = 1e-5; 

switch adaptivity_data.est_type
    case 0, estimator_type = 'none'; adaptivity_data.flag = 'functions';
    case 1, estimator_type = 'estandar'; adaptivity_data.flag = 'elements';
    case 2, estimator_type = 'nonestandar'; adaptivity_data.flag = 'elements';
    case 3, estimator_type = 'basis_functions'; adaptivity_data.flag = 'functions';
    case 4, estimator_type = 'basis_functions_2'; adaptivity_data.flag = 'functions';
end

% GRAFICOS
ver = 1;
plot_discrete_solution = ver;
plot_hierarchical_mesh = ver;
print_graphics = 0;
geps = 0;
pausas = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the hierarchical mesh and space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hmsh, hspace, geometry] = adaptivity_initialize (problem_data, method_data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Printing the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmpi (method_data.space_type, 'simplified'))
  space_type = 0;
else
  space_type = 1;
end

outputfile = sprintf('output_files/aigmdata_problem%d_degree%d_est%d_mark%2d_basis%d.txt',problem,method_data.degree(1),...
    adaptivity_data.est_type,adaptivity_data.mark_param*100,space_type);
file = fopen(outputfile,'w');
fprintf(file,'%% DOFs   elem   H1-err   EST     EST/H1-err    NBF    NoL    h_max \n');
fclose(file);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADAPTIVE LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s\n',outputfile);
% indices = [];
iter = 0;

if adaptivity_data.num_max_iter <= 0
    return
end

while 1
    
    iter = iter + 1;
    
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n',iter);
    fprintf('%s\n',outputfile);
    
    if ~hspace_check_partition_of_unity(hspace, hmsh)
        disp('ERROR: The partition-of-the-unity property does not hold.')
        return,
    end
    
    if (plot_hierarchical_mesh)
        hmsh_plot_cells (hmsh, 1);
%        plot_hmesh_param(hmsh, 1); % In figure(1)
        if print_graphics
            filename = sprintf('-problem%d-degree%d-est%d-iter%03d-',problem, method_data.degree(1),adaptivity_data.est_type,iter);
            print('-dpng', ['meshes/mesh' filename ])
            if geps
                print('-depsc2', ['meshes/mesh' filename])
            end
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SOLVE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    u = adaptivity_solve (hmsh, hspace, problem_data);
    
    if plot_discrete_solution
       plot_numerical_and_exact_solution(u, hmsh, hspace, problem_data.uex), 
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Error computation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [err_h1, err, err_h1s] = hspace_h1_error (hspace, hmsh, u, problem_data.uex, problem_data.graduex);
    fprintf('error_H1s = %g\n', err_h1s);
    
    clear err_h1_loc err_loc err_h1s_loc
    
    [~, ~, seminorma_h1s] = hspace_h1_error (hspace, hmsh, zeros (size(u)), problem_data.uex, problem_data.graduex);
    fprintf('seminorma_H1s = %g\n', seminorma_h1s);
    clear aaa
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ESTIMATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    disp('Computing error estimators:')
    est = estimate_laplace_param(u, hmsh, hspace, problem_data, adaptivity_data);
    gest = norm(est);
    tempo = toc;
    fprintf('EST: %f (%f seconds)\n', gest, tempo);
    
    hh = hmsh_get_element_size (hmsh);
    hh = hh(1);
    
    file = fopen(outputfile,'a');
    fprintf(file,'%5d  %5d %2.8f  %2.8f  %2.4f    %2d %2d  %.6f\n',...
        hspace.ndof, hmsh.nel, err_h1s, gest, gest/err_h1s, 0, hspace.nlevels, hh);
    fclose(file);
   
    
    if gest < adaptivity_data.tol
        disp('Success: The error estimation reached the desired tolerance')
        return,
    end
    
    if iter == adaptivity_data.num_max_iter
        disp('Warning: Maximum amount of iterations reached')
        return;
    end
    
    if hmsh.nlevels >= adaptivity_data.max_level
        disp('Warning: Maximum amount of levels reached')
        return;
    end
    
    if hspace.ndof > adaptivity_data.max_ndof
        disp('Warning: Maximum allowed DOFs achieved')
        return;
    end
    
    if hmsh.nel > adaptivity_data.max_nel
        disp('Warning: Maximum allowed amount of elements achieved')
        return;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MARK
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tic
    disp('Marking:')
    [marked, num_marked] = adaptivity_mark(est, hmsh, hspace, adaptivity_data);
    tempo = toc;
    nest = numel(est);
    switch adaptivity_data.flag
        case 'elements',
            fprintf('%d elements marked for refinement (%f seconds)\n', num_marked, tempo);
        case 'functions',
            fprintf('%d functions marked for refinement (%f seconds)\n', num_marked, tempo);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% REFINE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
    
    
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');
    if pausas
    pause
    end
end
