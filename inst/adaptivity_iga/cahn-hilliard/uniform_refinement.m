%--------------------------------------------------------------------------
% multiple refinements marking all the elements 
%--------------------------------------------------------------------------
function [hmsh, hspace] = uniform_refinement(hmsh, hspace, n_refinements, adaptivity_data)

  if (n_refinements >= 1)
    for ii = 1:n_refinements
      if (hmsh.nlevels >= adaptivity_data.max_level) % check max depth
        disp('Uniform refinement limited by max_level')
        break
      end

      marked = cell (hmsh.nlevels,1);
      marked{hmsh.nlevels} = hmsh.active{hmsh.nlevels};
      [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
    end
  end
end
