% Adaptive isogeometric method for d-dimensional problems (d = 1,2,3)

clear hmsh hspace u est
% close all
add_my_paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1
    problem = 6; 
    est_type = 3;
    mark_param = .5;
    p = 3; % polynomial degree
    flag_whole_basis = 0;
end

ver = 0;
plot_discrete_solution = ver;
plot_hierarchical_mesh = ver;
print_graphics = 0;
geps = 0;
pausas = 0;

mark_strategy = 'MS';
max_level = 10;
max_ndof = 5000;
num_max_iter = 15;
max_nel = 5000;
% dim = 2; % Esto segun el problema
initial_num_el = 2;

switch est_type
    case 0, estimator_type = 'none'; flag = 'functions';
    case 1, estimator_type = 'estandar'; flag = 'elements';
    case 2, estimator_type = 'nonestandar'; flag = 'elements';
    case 3, estimator_type = 'basis_functions'; flag = 'functions';
    case 4, estimator_type = 'basis_functions_2'; flag = 'functions';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[problem_data, exact_solution] = load_problem_data(problem);

dim = problem_data.dim;

uex = exact_solution.uex;
graduex = exact_solution.graduex;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the hierarchical mesh and space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tp_msh, tp_space] = get_initial_msh_and_space(dim, p, initial_num_el, problem_data.geo_name);
[hmsh, hspace] = tp2hier (tp_msh, tp_space, problem_data.geo_name, flag_whole_basis);

clear tp_msh tp_space


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Printing the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outputfile = sprintf('output_files/aigmdata_problem%d_degree%d_est%d_mark%2d_basis%d.txt',problem,p,est_type,mark_param*100,flag_whole_basis);
file = fopen(outputfile,'w');
fprintf(file,'%% DOFs   elem   H1-err   EST     EST/H1-err    NBF    NoL    h_max \n');
fclose(file);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADAPTIVE LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s\n',outputfile);
% indices = [];
iter = 0;

if num_max_iter <= 0
    return
end

while 1
    
    iter = iter + 1;
    
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n',iter);
    fprintf('%s\n',outputfile);
    
    if ~check_partition_of_the_unity(hmsh, hspace)
        disp('ERROR: The partition-of-the-unity property does not hold.')
        return,
    end
    
    if (plot_hierarchical_mesh && hmsh.ndim > 1)
        plot_msh(hmsh, iter); % In figure(1)
        if print_graphics
            filename = sprintf('-problem%d-degree%d-est%d-iter%03d-',problem, p,est_type,iter);
            print('-dpng', ['meshes/mesh' filename ])
            if geps
                print('-depsc2', ['meshes/mesh' filename])
            end
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SOLVE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    u = assemble_and_solve(hmsh, hspace, problem_data);
    
    if plot_discrete_solution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot numerical solution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        switch dim
            case 1, %x = hmsh.quad_nodes(:);
                [Z, x] = hspline_eval(u, hmsh, hspace, 4);
                figure(2)
                subplot (1,2,1)
                plot (x(:), Z(:),'*')
                title ('Numerical solution'), axis tight
                subplot (1,2,2)
                x= linspace(0,1);
                u_ex_values = uex(x);
                plot (x, u_ex_values)
                title ('Exact solution'), axis tight
                
            case 2,
                plot_numerical_and_exact_solution(u, hmsh, hspace, uex) 
        end
        clear Z ZZ x y X Y 
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Error computation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [err_h1, err, err_h1s, err_h1_loc, err_loc, err_h1s_loc] = compute_H1_error(u, uex, graduex, hmsh, hspace);
    fprintf('error_H1s = %g\n', err_h1s);
    
    clear err_h1_loc err_loc err_h1s_loc
    
    [aaa, aaa, seminorma_h1s, aaa, aaa, aaa] = compute_H1_error(zeros(hspace.ndof,1), uex, graduex, hmsh, hspace);
    fprintf('seminorma_H1s = %g\n', seminorma_h1s);
    clear aaa
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ESTIMATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    disp('Computing error estimators:')
    est = estimate(u, hmsh, hspace, problem_data, estimator_type);
    gest = norm(est);
    tempo = toc;
    fprintf('EST: %f (%f seconds)\n', gest, tempo);
    
    hh = get_meshsize(hmsh);
    hh = hh(1);
    
    file = fopen(outputfile,'a');
    fprintf(file,'%5d  %5d %2.8f  %2.8f  %2.4f    %2d %2d  %.6f\n',...
        hspace.ndof, hmsh.nel, err_h1s, gest, gest/err_h1s, 0, hspace.nlevels, hh);
    fclose(file);
   
    
    if gest < 1e-8
        disp('Success: The error estimation reached the desired tolerance')
        return,
    end
    
    if iter == num_max_iter
        disp('Warning: Maximum amount of iterations reached')
        return;
    end
    
    if hmsh.nlevels >= max_level
        disp('Warning: Maximum amount of levels reached')
        return;
    end
    
    if hspace.ndof > max_ndof
        disp('Warning: Maximum allowed DOFs achieved')
        return;
    end
    
    if hmsh.nel > max_nel
        disp('Warning: Maximum allowed amount of elements achieved')
        return;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MARK
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tic
    disp('Marking:')
    [marked, num_marked] = mark(est, hmsh, hspace, mark_strategy, mark_param,flag);
    tempo = toc;
    nest = numel(est);
    switch flag
        case 'elements',
            fprintf('%d elements marked for refinement (%f seconds)\n', num_marked, tempo);
        case 'functions',
            fprintf('%d functions marked for refinement (%f seconds)\n', num_marked, tempo);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% REFINE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [hmsh, hspace] = refine (hmsh, hspace, marked, flag);
    
    
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');
    if pausas
    pause
    end
end
