function  adaptivity_solve_laplace(problem_data, method_data, adaptivity_data)


[hmsh, hspace] = init_hierarchical_mesh_and_space(problem_data,method_data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the hierarchical mesh and space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tp_msh, tp_space] = get_initial_msh_and_space(dim, p, initial_num_el, problem_data.geo_name);
[hmsh, hspace] = tp2hier (tp_msh, tp_space, problem_data.geo_name, flag_whole_basis);

clear tp_msh tp_space


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
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Error computation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [err_h1, err, err_h1s, err_h1_loc, err_loc, err_h1s_loc] = compute_H1_error(u, uex, graduex, hmsh, hspace);
    fprintf('error_H1s = %g\n', err_h1s);
    
    clear err_h1_loc err_loc err_h1s_loc
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ESTIMATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     tic
%     disp('Computing error estimators:')
    est = estimate(u, hmsh, hspace, problem_data, estimator_type);
    gest = norm(est);
%     tempo = toc;
%     fprintf('EST: %f (%f seconds)\n', gest, tempo);
    
%     hh = get_meshsize(hmsh);
%     hh = hh(1);
%     
%    
    
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
    
   % tic
   % disp('Marking:')
    [marked, num_marked] = mark(est, hmsh, hspace, mark_strategy, mark_param,flag);
%     tempo = toc;
%     nest = numel(est);
%     switch flag
%         case 'elements',
%             fprintf('%d elements marked for refinement (%f seconds)\n', num_marked, tempo);
%         case 'functions',
%             fprintf('%d functions marked for refinement (%f seconds)\n', num_marked, tempo);
%     end
%     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% REFINE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [hmsh, hspace] = refine (hmsh, hspace, marked, flag);
    
    
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');
    
end
