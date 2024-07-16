function est = adaptivity_estimate_cahn_hilliard (u, hmsh, hspace, adaptivity_data)

  % gradient-based estimator ------------------------------------------------
  if strcmp(adaptivity_data.estimator_type, 'gradient') == 1
    normalization = 1;

    % compute the average norm of the gradient 
    integral = integrals_elem (hspace, hmsh, u, 'norm_grad');        
    est = integral /normalization;

  % field-based estimator ---------------------------------------------------
  elseif strcmp(adaptivity_data.estimator_type, 'field') == 1

    % compute the average value of the field       
    integral = integrals_elem (hspace, hmsh, u, 'field');            
    est = 1-(abs(integral));
  end
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

      if strcmp(type, 'field') == 1
        sp_level = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
      elseif strcmp(type, 'norm_grad') == 1
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

  if (strcmp(type, 'field') == 1)
    eu = sp_eval_msh (u, sp, msh,'value');
    valu = reshape (eu, msh.nqn, msh.nel);
    integral = (sum (valu .* w)) ./ elem_size;

  elseif (strcmp(type, 'norm_grad') == 1)
    eu = sp_eval_msh (u, sp, msh,'gradient');
    grad_valu = reshape (eu, sp.ncomp, msh.rdim, msh.nqn, msh.nel);
    integral = sum (reshape ( sqrt (sum ((grad_valu).^2, 2) ), [msh.nqn, msh.nel]) .* w);
    integral =  integral ./ elem_size;
  end
end
