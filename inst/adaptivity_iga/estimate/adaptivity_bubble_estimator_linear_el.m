% ADAPTIVITY_BUBBLE_ESTIMATOR_LINEAR_EL: computation of a posteriori error indicators for linear elasticity problem, using bubble functions
%
% USAGE:
%
%   est = adaptivity_bubble_estimator_linear_el (u, hmsh, hspace, problem_data)
%
% INPUT:
%
%   u:            degrees of freedom
%   hmsh:         object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace:       object representing the space of hierarchical splines (see hierarchical_space)
%   problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - lambda_lame:   function handle of the first Lame' parameter
%    - mu_lame:       function handle of the second Lame' parameter
%    - f:             function handle of the source term
%   adaptivity_data: structure with data for adaptivity. For this function, it must contain the field:
%    - C0_est:         multiplicative constant for the error indicators
%
%
% OUTPUT:
%
%   est: computed a posteriori error indicators 
%
% WARNING: the current version of the code does not modify the bubble space
%           depending on the boundary conditions
%
%   For more information on the bubble estimators, see Coradello, Antolin and Buffa, CMAME, 2020.
%         
%
% Copyright (C) 2018-2023 Luca Coradello, Rafael Vazquez
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

function estimator = adaptivity_bubble_estimator_linear_el (u, hmsh, hspace, problem_data, adaptivity_data)
  
  if (nargin < 5 || ~isfield (adaptivity_data, 'C0_est'))
    C0_est = 1;
  else
    C0_est = adaptivity_data.C0_est;
  end
  lambda = problem_data.lambda_lame;
  mu = problem_data.mu_lame;

  last_dof = cumsum (hspace.ndof_per_level);

  if (isa (hspace.space_of_level(1), 'sp_scalar'))
    error ('Not implemented yet')
    deg = hspace.space_of_level(1).degree;
    space_bubble = space_bubble_function_bilaplacian (hmsh, deg);
    nqn = hmsh.mesh_of_level(1).nqn;
  elseif (isa (hspace.space_of_level(1), 'sp_multipatch_C1'))
    deg = hspace.space_of_level(1).sp_patch{1}.degree;    
    space_bubble = space_bubble_function_bilaplacian_mp (hmsh, deg);
    nqn = hmsh.mesh_of_level(1).msh_patch{1}.nqn;
  else
    error ('Estimator only implemented for multipatch C1 scalar space objects')
  end  
  
  estimator = zeros(hmsh.nel, 1);
  
  shifting_vector = cumsum([0 hmsh.nel_per_level]);
  
  for ilev = 1:hmsh.nlevels
    if (hmsh.nel_per_level(ilev) > 0)
      x = cell (hmsh.rdim,1);
      for idim = 1:hmsh.rdim
        x{idim} = reshape (hmsh.msh_lev{ilev}.geo_map(idim,:,:), nqn, hmsh.nel_per_level(ilev));
      end
      spu_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true, 'gradient', true);
      spu_lev = change_connectivity_localized_Csub (spu_lev, hspace, ilev);
      spv_lev = space_bubble.space_of_level(ilev);
      spu_lev = sp_scalar2vector (spu_lev, hmsh.msh_lev{ilev}, 'value', false, 'gradient', true, 'divergence', true);
      spv_lev = sp_scalar2vector (spv_lev, hmsh.msh_lev{ilev}, 'value', true, 'gradient', true, 'divergence', true);
      
      dofs_to_lev = [];
      for icomp = 1:hmsh.rdim
        dofs_to_lev = union (dofs_to_lev, (icomp-1)*hspace.ndof + (1:last_dof(ilev)));
      end

      u_lev = reshape (u(dofs_to_lev), last_dof(ilev), hmsh.rdim);
      solution_of_level = reshape (hspace.Csub{ilev}*u_lev, [], 1);
      K_lev = op_su_ev (spu_lev, spv_lev, hmsh.msh_lev{ilev}, lambda(x{:}), mu(x{:}));
      b_lev = op_f_v (spv_lev, hmsh.msh_lev{ilev}, problem_data.f(x{:}));           
      residual_of_level = b_lev - K_lev * solution_of_level;

      K_err_lev = op_su_ev (spv_lev, spv_lev, hmsh.msh_lev{ilev}, lambda(x{:}), mu(x{:}));
      error_of_level = K_err_lev \ residual_of_level;

      conn = arrayfun (@(x) spv_lev.connectivity(1:spv_lev.nsh(x), x), 1:hmsh.msh_lev{ilev}.nel, 'UniformOutput', false);
      err_elem = cellfun(@(ind) full (error_of_level(ind).' * K_err_lev(ind,ind) * error_of_level(ind)), conn);
      err_elem = sqrt (err_elem);
      
      estimator((shifting_vector(ilev)+1):shifting_vector(ilev+1)) = err_elem;
    end
  end

  estimator = estimator * C0_est;

end