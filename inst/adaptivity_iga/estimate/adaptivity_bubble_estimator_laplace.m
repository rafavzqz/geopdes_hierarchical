% ADAPTIVITY_BUBBLE_ESTIMATOR_LAPLACE: computation of a posteriori error indicators for Laplacian problem, using bubble functions
%
% USAGE:
%
%   est = adaptivity_bubble_estimator_laplace (u, hmsh, hspace, problem_data)
%
% INPUT:
%
%   u:            degrees of freedom
%   hmsh:         object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace:       object representing the space of hierarchical splines (see hierarchical_space)
%   problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - c_diff:        diffusion coefficient (epsilon in the equation), assumed to be a smooth (C^1) function
%    - f:             function handle of the source term
%   adaptivity_data: structure with data for adaptivity. For this function, it must contain the field:
%    - C0_est:         multiplicative constant for the error indicators
%
%
% OUTPUT:
%
%   est: computed a posteriori error indicators 
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

function estimator = adaptivity_bubble_estimator_laplace (u, hmsh, hspace, problem_data, adaptivity_data)
  
  if(~isfield(problem_data, 'c_diff'))  
    problem_data.c_diff  = @(x, y) ones(size(x));
  end
  if (nargin < 5 || ~isfield (adaptivity_data, 'C0_est'))
    C0_est = 1;
  else
    C0_est = adaptivity_data.C0_est;
  end
  
  last_dof = cumsum (hspace.ndof_per_level);

  if (isa (hspace.space_of_level(1), 'sp_scalar'))
    deg = hspace.space_of_level(1).degree;
    space_bubble = space_bubble_function_laplacian (hmsh, deg);
    nqn = hmsh.mesh_of_level(1).nqn;
  elseif (isa (hspace.space_of_level(1), 'sp_multipatch'))
    deg = hspace.space_of_level(1).sp_patch{1}.degree;    
    space_bubble = space_bubble_function_laplacian_mp (hmsh, deg);
    nqn = hmsh.mesh_of_level(1).msh_patch{1}.nqn;
  end  
  
  estimator = zeros(1,hmsh.nel);
  
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
      
      solution_of_level = hspace.Csub{ilev}*u(1:last_dof(ilev));
      K_lev = op_gradu_gradv (spu_lev, spv_lev, hmsh.msh_lev{ilev}, problem_data.c_diff (x{:}));      
      b_lev = op_f_v (spv_lev, hmsh.msh_lev{ilev}, problem_data.f(x{:}));           
      residual_of_level = b_lev - K_lev * solution_of_level;

      K_err_lev = op_gradu_gradv (spv_lev, spv_lev, hmsh.msh_lev{ilev}, problem_data.c_diff (x{:}));
      error_of_level = K_err_lev \ residual_of_level;
      
% Compute the estimator from the matrix, avoiding a loop on the elements
      conn = arrayfun (@(x) spv_lev.connectivity(1:spv_lev.nsh(x), x), 1:hmsh.msh_lev{ilev}.nel, 'UniformOutput', false);
      err_elem = cellfun(@(ind) full(error_of_level(ind).' * K_err_lev(ind,ind) * error_of_level(ind)), conn);
      err_elem = sqrt (err_elem);
      
      estimator((shifting_vector(ilev)+1):shifting_vector(ilev+1)) = err_elem;
    end
  end
  
  estimator = estimator * C0_est;

end