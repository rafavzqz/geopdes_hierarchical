% TEST

% PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));
problem_data.grad_c_diff = @(x, y) cat (1, ...
    reshape (zeros(size(x)), [1, size(x)]), ...
    reshape (zeros(size(x)), [1, size(x)]));

% Source and boundary terms
% Choose 1.5 < Cx, Cy < 2.5. Thus, the solution belongs H2\H3
Cx = 1.6;
Cy = 2.4;
problem_data.f = @(x,y) Cx*((Cx+1)*x.^(Cx-1)-(Cx-1)*x.^(Cx-2)).*y.^Cy.*(1-y)+...
    Cy*((Cy+1)*y.^(Cy-1)-(Cy-1)*y.^(Cy-2)).*x.^Cx.*(1-x);
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) zeros(size(x));

% Exact solution (optional)
problem_data.uex = @(x,y) x.^Cx.*(1-x).*y.^Cy.*(1-y);
problem_data.graduex = @(x,y) cat (1, ...
    reshape ((Cx*x.^(Cx-1).*(1-x)-x.^Cx).*y.^Cy.*(1-y), [1, size(x)]), ...
    reshape ((Cy*y.^(Cy-1).*(1-y)-y.^Cy).*x.^Cx.*(1-x), [1, size(x)]));


% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = [3 3];        % Degree of the splines
method_data.regularity  = [2 2];        % Regularity of the splines
method_data.nsub_coarse = [2 2];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];        % Number of subdivisions for each refinement
method_data.nquad       = [4 4];        % Points for the Gaussian quadrature rule
%method_data.space_type  = 'simplified'; % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 0;            % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
%adaptivity_data.flag = 'elements';
%adaptivity_data.flag = 'functions';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = .5;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level = 15;
adaptivity_data.max_ndof = 10000;
%adaptivity_data.num_max_iter = 16;
adaptivity_data.max_nel = 10000;
adaptivity_data.tol = 1e-5;

for test_ind = 1:3
    switch test_ind
        case 1, adaptivity_data.flag = 'elements';
            method_data.space_type  = 'simplified';
            adaptivity_data.num_max_iter = 13;
        case 2, adaptivity_data.flag = 'elements';
            method_data.space_type  = 'standard';
            adaptivity_data.num_max_iter = 13;
        case 3, adaptivity_data.flag = 'functions';
            method_data.space_type  = 'simplified';
            adaptivity_data.num_max_iter = 14;
    end
    
    outputfile = sprintf('test_result_%d.txt',test_ind);
    file = fopen(outputfile,'w');
    fprintf(file,'%% DOFs   elem   H1-err   EST     EST/H1-err\n');
    fclose(file);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Definition of the hierarchical mesh and space
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('%s\n',outputfile);
    iter = 0;
    [hmsh, hspace] = adaptivity_initialize_laplace (problem_data, method_data);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ADAPTIVE LOOP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    while 1
        
        iter = iter + 1;
        
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n',iter);
        fprintf('%s\n',outputfile);
        
        if ~hspace_check_partition_of_unity(hspace, hmsh)
            disp('ERROR: The partition-of-the-unity property does not hold.')
            return,
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% SOLVE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        u = adaptivity_solve_laplace (hmsh, hspace, problem_data);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Error computation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [~,~,err_h1s] = sp_h1_error (hspace, hmsh, u, problem_data.uex, problem_data.graduex);
        fprintf('error_H1s = %g\n', err_h1s);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ESTIMATE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic
        disp('Computing error estimators:')
        est = adaptivity_estimate_laplace(u, hmsh, hspace, problem_data, adaptivity_data);
        gest = norm(est);
        tempo = toc;
        fprintf('EST: %f (%f seconds)\n', gest, tempo);
        
        file = fopen(outputfile,'a');
        fprintf(file,'%5d  %5d %2.8f  %2.8f  %2.4f\n',...
            hspace.ndof, hmsh.nel, err_h1s, gest, gest/err_h1s);
        fclose(file);
        
        if iter == adaptivity_data.num_max_iter
            disp('Warning: Maximum amount of iterations reached')
            break;
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
        
    end
    
    
end