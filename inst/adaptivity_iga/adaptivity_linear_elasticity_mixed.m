% ADAPTIVITY_linear_elasticity: solve the linear elasticity problem with an adaptive isogeometric method based on hierarchical splines.
%

function [geometry, hmsh, hspace, u, hspace_press, press, solution_data] = adaptivity_linear_elasticity_mixed_adma2 (problem_data, method_data, adaptivity_data, plot_data, savename_txt)


if (nargin == 3)
  plot_data = struct ('print_info', true, 'plot_hmesh', false, 'plot_discrete_sol', false);
end
if (~isfield (plot_data, 'print_info'))
  plot_data.print_info = true;
end
if (~isfield (plot_data, 'plot_hmesh'))
  plot_data.plot_hmesh = false;
end
if (~isfield (plot_data, 'plot_discrete_sol'))
  plot_data.plot_discrete_sol = false;
end

% Initialization of some auxiliary variables
if (plot_data.plot_hmesh)
  fig_mesh = figure;
end
if (plot_data.plot_discrete_sol)
  fig_sol = figure;
end
nel = zeros (1, adaptivity_data.num_max_iter); ndof = nel; gest = nel+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the hierarchical mesh and space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hmsh, hspace, hspace_press, geometry] = adaptivity_initialize_TH (problem_data, method_data);

global outputfile file;
%outputfile = sprintf('exm02bD2adm0.txt');
%outputfile = sprintf('exm02bD2adm2.txt');
%outputfile = sprintf('exm02bD2adm3.txt');
outputfile = sprintf(savename_txt);
file = fopen(outputfile,'w');
fprintf(file,'%%--------------------------------------------------------------------------------------------------------------------%%\n');
fprintf(file,'%% iter l(space)  DOFs elem    H1-err        EST     EI(EST/H1-err)   EI_min       EI_max\n      sparsity       bandwidth      conditioning');
fprintf(file,'%%--------------------------------------------------------------------------------------------------------------------%%\n');
fclose(file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADAPTIVE LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter = 0;
while (iter < adaptivity_data.num_max_iter)
  iter = iter + 1;
  
  if (plot_data.print_info)
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',iter);
  end
  if (~hspace_check_partition_of_unity (hspace, hmsh))
    disp('ERROR: The partition-of-the-unity property does not hold.')
    solution_data.flag = -1; break
  end
  if (~hspace_check_partition_of_unity (hspace_press, hmsh))
    disp('ERROR: The partition-of-the-unity property does not hold.')
    solution_data.flag = -1; break
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% SOLVE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  curr_lev(iter)=hspace.nlevels;
  curr_dof(iter)=hspace.ndof;
  curr_el(iter)=hmsh.nel;  
  if (plot_data.print_info)
    disp('SOLVE:')
    fprintf('Number of levels: %d \n', hspace.nlevels)
    fprintf('Number of levels (mesh): %d \n', hmsh.nlevels)
    fprintf('Number of elements: %d. Total DOFs (space 1): %d. Total DOFs (space 2): %d \n', hmsh.nel, hspace.ndof, hspace_press.ndof);
  end
  [u,press,condnum(iter),sparsity(iter,:),bandwidth(iter)] = adaptivity_solve_linear_elasticity_mixed_plus (hmsh, hspace, hspace_press, problem_data);  %add outputs for condition number, sparsity e bandwidth?
  nel(iter) = hmsh.nel; ndof(iter) = hspace.ndof;  ndof2(iter) = hspace_press.ndof;
  
  displ = sp_eval (u, hspace, geometry, {1-eps, 1-eps});
  displace(iter) = displ(2);
  
  if (plot_data.plot_hmesh)
    hmsh_plot_cells (hmsh, 10, fig_mesh);
  end
  if (plot_data.plot_discrete_sol)
    figure(fig_sol)  
    pts = {linspace(0, 1, 21), linspace(0, 1, 21)};  
    [eu, F] = sp_eval (u, hspace, geometry, pts);
    [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
    subplot (1,2,1)
    quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
    title ('Numerical solution'), axis equal tight
    subplot (1,2,2)
    eu2 = problem_data.uex (X, Y);
    quiver (X, Y, squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
    title ('Exact solution'), axis equal tight  
  end
  if (plot_data.plot_hmesh || plot_data.plot_discrete_sol)
%     disp('Paused. Type "dbcont" to continue')
%     keyboard
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% ESTIMATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp('ESTIMATE:')
   if (plot_data.print_info); disp('ESTIMATE:'); end
     est = adaptivity_estimate_linear_el_mixed(u, press, hmsh, hspace, hspace_press, problem_data, adaptivity_data);
     est_elem{iter} = est.';
     gest(iter) = norm (est);
   if (plot_data.print_info); fprintf('Computed estimate: %e \n', gest(iter)); end
   if (isfield (problem_data, 'uex'))
     [~, errl2(iter), err_h1s(iter),~,errl2_elem{iter},err_h1s_elem{iter}] = sp_h1_error (hspace, hmsh, u, problem_data.uex, problem_data.graduex);
     if (plot_data.print_info); fprintf('Error in H1 seminorm = %e\n', err_h1s(iter)); end
     %efficiency indicator
     if strcmp(adaptivity_data.flag,'elements') & (size(est_elem{iter})==size(err_h1s_elem{iter}))
       efficiency_elem{iter}=est_elem{iter}./err_h1s_elem{iter};
     else
       efficiency_elem{iter}=zeros(hmsh.nel,1);  
     end
     efficiency(iter)=gest(iter)/err_h1s(iter);
     if (plot_data.print_info)
       fprintf('Efficiency = %g\n', efficiency(iter));
       fprintf('Min efficiency (elements) = %g\n', min(efficiency_elem{iter}));
       fprintf('Max efficiency (elements) = %g\n', max(efficiency_elem{iter}));
     end
     
     file = fopen(outputfile,'a');
     fprintf(file,'%5d & %5d & %5d & %5d & %5e & %5e & %e & %e & %e \\ \n',...
     iter, hspace.nlevels, hspace.ndof, hmsh.nel, err_h1s(iter), gest(iter),...
     efficiency(iter), min(efficiency_elem{iter}), max(efficiency_elem{iter}), sparsity(iter), bandwidth(iter), condnum(iter));
     fclose(file);   
   end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %% STOPPING CRITERIA
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (gest < adaptivity_data.tol)
    disp('Success: The error estimation reached the desired tolerance')
    break
  end
    
  if (iter == adaptivity_data.num_max_iter)
    disp('Warning: Maximum amount of iterations reached')
    break
  end
    
  if (hmsh.nlevels >= adaptivity_data.max_level)
    disp('Warning: Maximum amount of levels reached')
    break
  end
    
  if (hspace.ndof > adaptivity_data.max_ndof)
    disp('Warning: Maximum allowed DOFs achieved')
    break
  end
    
  if (hmsh.nel > adaptivity_data.max_nel)
    disp('Warning: Maximum allowed amount of elements achieved')
    break
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% MARK
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (plot_data.print_info); disp('MARK:'); end
  [marked, num_marked] = adaptivity_mark (est, hmsh, hspace, adaptivity_data);

%   for i=1:hmsh.nlevels-1
%       marked{i}=[];
%   end
%   marked{hmsh.nlevels}=1:hmsh.nel_per_level(hmsh.nlevels);
%   num_marked=hmsh.nel;
  if (plot_data.print_info); 
    fprintf('%d %s marked for refinement \n', num_marked, adaptivity_data.flag);
    disp('REFINE:')
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% REFINE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [hmsh, hspace, hspace_press] = adaptivity_refine_TH_adm (hmsh, hspace, hspace_press, marked, adaptivity_data);
end

solution_data.iter = iter;
solution_data.gest = gest(1:iter);
solution_data.ndof = ndof(1:iter);
solution_data.ndof2 = ndof2(1:iter);
solution_data.nel  = nel(1:iter);
solution_data.displace = displace(1:iter);
if (exist ('condnum', 'var')) 
  solution_data.condnum  = condnum(1:iter);
end
if (exist ('sparsity', 'var')) 
  solution_data.sparsity  = sparsity(1:iter,:);
end
if (exist ('bandwidth', 'var')) 
  solution_data.bandwidth = bandwidth(1:iter);
end
if (exist ('err_h1s', 'var')) 
  solution_data.err_h1s = err_h1s(1:iter);
  solution_data.err_h1s_elem=err_h1s_elem;
end
if (exist('gest','var'))
  solution_data.est=gest(1:iter);
  solution_data.est_elem=est_elem;
end
if (exist('efficiency','var'))
  solution_data.efficiency=efficiency(1:iter);
  solution_data.efficiency_elem=efficiency_elem;
end

end
