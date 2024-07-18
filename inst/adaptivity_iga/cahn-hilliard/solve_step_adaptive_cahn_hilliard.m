% SOLVE_STEP_ADAPTIVE_CAHN_HILLIARD: perform one step of the generalized alpha method
%  for the solution of the Cahn-Hilliard equation, and adapts the mesh,
%  doing both refinement and coarsening.
%  It is called from adaptive_cahn_hilliard or adaptive_cahn_hilliard_mp_C1
%
% INPUT:
%
%  u_n:          field at the previous time step.
%  u_dotn:       time derivative at the previous time step.
%  hspace:       space object (see hierarchical_space or hierarchical_space_mp_C1)
%  hmsh:         mesh object (see hierarchical_mesh or hierarchical_mesh_mp)
%  dt:           time step size
%  a_m:          parameter for the generalized alpha method
%  a_f:          parameter for the generalized alpha method
%  gamma:        parameter for the generalized alpha method
%  Cpen:         penalization parameter for Nitsche's method
%  problem_data: a structure with data of the problem.
%  adaptivity_data: a structure with data for adaptivity (in space).
%  old_space:    space from previous iterations. If not changed, some matrices are not recomputed.
%  nmnn_sides:   sides where to impose the Neumann condition.
%
% OUTPUT:
%
%  u_n1:      field at the new time step.
%  u_dotn1:   time derivative at the new time step.
%  hspace:    space object, it can be refined or coarsened (see hierarhical_space or hierarchical_space_mp_C1)
%  hmsh:      mesh object, it can be refined or coarsened (see hierarchical_mesh or hierarchical_mesh_mp)
%  est:       refinement indicator
%  old_space: space of the previous iteration, updated
%
% Copyright (C) 2023, 2024 Michele Torre, Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [u_n1, udot_n1, hspace, hmsh, est, old_space] = solve_step_adaptive_cahn_hilliard...
  (u_n, udot_n, hspace, hmsh, dt, a_m, a_f, gamma, Cpen, ...
   problem_data, adaptivity_data, old_space, nmnn_sides)

  lambda = problem_data.lambda;
  mu = problem_data.mu;
  dmu = problem_data.dmu;
  iter = 0;
  while (1)
    iter = iter + 1;
    disp(strcat('%%%%%%%%%%%%%%%%% Adaptivity iteration ',num2str(iter),' %%%%%%%%%%%%%%%%%'));
        
    %------------------------------------------------------------------
    % solve
    [u_n1, udot_n1, old_space] = generalized_alpha_step_cahn_hilliard...
                                  (u_n, udot_n, dt, a_m, a_f, gamma, lambda, mu, dmu, ...
                                   Cpen, hspace, hmsh, old_space, nmnn_sides);

    %------------------------------------------------------------------
    %estimate
% %%%%%%%%%%%%%%%%%%%%%%%%%% TODO: write the help of the function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    est = adaptivity_estimate_cahn_hilliard (u_n1, hmsh, hspace, adaptivity_data);

    %------------------------------------------------------------------
    % stopping criteria
    if (iter == adaptivity_data.num_max_iter)
      disp('Warning: reached the maximum number of iterations')
    elseif (hmsh.nlevels > adaptivity_data.max_level)
      disp(strcat('number of levels =',num2str(hmsh.nlevels))) 
      disp('Warning: reached the maximum number of levels')
      break
    end
        
    if (hmsh.nlevels == adaptivity_data.max_level && hmsh.nel == hmsh.nel_per_level(end))
      disp('Mesh completely refined at the maximum level')
      break
    end

    %------------------------------------------------------------------
    % mark
% %%%%%%%%%%%%%%%%%%%%%%%%%% TODO: write the help of the function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    [marked, num_marked_ref] = adaptivity_mark_cahn_hilliard (est, hmsh, hspace, adaptivity_data);
    % limit the maximum refinement depth
    if (hmsh.nlevels == adaptivity_data.max_level)
      num_deleted = numel(marked{hmsh.nlevels});
      marked{hmsh.nlevels} = [];
      num_marked_ref = num_marked_ref - num_deleted ;
    end

    % stopping criterion
    if (num_marked_ref == 0)
      disp('No element is refined')
      break  
    end

    % refine
    [hmsh, hspace, Cref] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
    old_space = struct ('modified', true, 'space', [], 'mesh', [], 'mass_mat', [], ...
                        'lapl_mat', [], 'bnd_mat', [], 'Pen', [], 'pen_rhs', []);

    % recompute control variables
    u_n = Cref * u_n;       
    udot_n = Cref * udot_n;
   
  end % end loop adaptivity
end
