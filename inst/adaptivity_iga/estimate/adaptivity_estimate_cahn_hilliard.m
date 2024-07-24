% ADAPTIVITY_ESTIMATE_CAHN_HILLIARD: Computation of an indication for Cahn-Hilliard equation, using globally smooth (C^1) hierarchical spaces.
%
% USAGE:
%
%   est = adaptivity_estimate_cahn_hilliard (u, hmsh, hspace, adaptivity_data)
%
% INPUT:
%
%   u:       degrees of freedom
%   hmsh:    object representing the hierarchical mesh (see hierarchical_mesh or hierarchical_mesh_mp)
%   hspace:  object representing the space of hierarchical splines (see hierarchical_space or hierarchical_space_mp_C1)
%   adaptivity_data: a structure with the data for the adaptivity method. In particular, it contains the fields:
%    - estimator_type: either 'field' or 'gradient'
%
%
% OUTPUT:
%
%   est: computed a posteriori error indicator
%
%  For the details, see the reference paper
%   C. Bracco, C. Giannelli, A. Reali, M. Torre, R. Vazquez, CMAME 417 (2023), 116355. 
%
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

function [est, mark_param_coarsening] = adaptivity_estimate_cahn_hilliard (u, hmsh, hspace, adaptivity_data)

  % gradient-based estimator ------------------------------------------------
  if (strcmpi(adaptivity_data.estimator_type, 'gradient'))

    % compute the average value of the gradient
    integral = integrals_elem (hspace, hmsh, u, 'norm_grad');        
    est = integral ;

  % field-based estimator ---------------------------------------------------
  elseif (strcmpi(adaptivity_data.estimator_type, 'field'))

    % compute the average value of the field
    integral = integrals_elem (hspace, hmsh, u, 'field');            
    est = 1 - (abs(integral));
  end

tmp = find(est < adaptivity_data.mark_param);
if numel(tmp)>0
    est(tmp) = 0;
end
est = est * adaptivity_data.mark_param;
mark_param_coarsening = numel(tmp)/numel(est);
end


%--------------------------------------------------------------------------
% functions
%--------------------------------------------------------------------------
function integral = integrals_elem (hspace, hmsh, u, type)

  if (numel(u) ~= hspace.ndof)
    error ('Wrong size of the vector of degrees of freedom')
  end

  integral = zeros (1, hmsh.nel);

  first_elem = cumsum ([0 hmsh.nel_per_level]) + 1;
  last_elem = cumsum ([hmsh.nel_per_level]);
  last_dof = cumsum (hspace.ndof_per_level);
  for ilev = 1:hmsh.nlevels
    if (hmsh.nel_per_level(ilev) > 0)
      msh_level = hmsh.msh_lev{ilev};

      if (strcmpi(type, 'field'))
        sp_level = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
      elseif (strcmpi(type, 'norm_grad'))
        sp_level = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', false, 'gradient', true);
      end

      sp_level = change_connectivity_localized_Csub (sp_level, hspace, ilev);
      int_elem = integrals_elem_level (sp_level, msh_level, hspace.Csub{ilev}*u(1:last_dof(ilev)), type);
    
      integral(:,first_elem(ilev):last_elem(ilev)) = int_elem;
    end
  end

end

%--------------------------------------------------------------------------
function integral = integrals_elem_level (sp, msh, u, type)

  w = msh.quad_weights .* msh.jacdet;
  elem_size = sum (ones(size(w)) .* w);

  if (strcmpi(type, 'field'))
    eu = sp_eval_msh (u, sp, msh, 'value');
    valu = reshape (eu, msh.nqn, msh.nel);
    integral = (sum (valu .* w)) ./ elem_size;

  elseif (strcmpi(type, 'norm_grad'))
    eu = sp_eval_msh (u, sp, msh, 'gradient');
    grad_valu = reshape (eu, sp.ncomp, msh.rdim, msh.nqn, msh.nel);
    integral = sum (reshape ( sqrt (sum ((grad_valu).^2, 2) ), [msh.nqn, msh.nel]) .* w);
    integral =  integral ./ elem_size;
  end
end
