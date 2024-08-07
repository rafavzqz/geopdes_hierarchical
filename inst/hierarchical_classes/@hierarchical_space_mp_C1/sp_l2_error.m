% SP_L2_ERROR: Evaluate the error in L^2 norm, for hierarchical splines.
%
%   [errl2, errl2_elem] = sp_l2_error (hspace, hmsh, u, uex)
%
% INPUT:
%
%    hspace:  object defining the discrete space (see hierarchical_space_mp_C1)
%    hmsh:    object representing the hierarchical mesh (see hierarchical_mesh_mp)
%    u:       vector of dof weights
%    uex:     function handle to evaluate the exact solution
%
% OUTPUT:
%
%     errl2:       error in L^2 norm
%     errl2_elem:  error in L^2 norm, for each single element
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2022-2023 Rafael Vazquez
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

function [errl2, errl2_elem] = sp_l2_error (hspace, hmsh, u, uex)

errl2 = 0;
errl2_elem = zeros (1, hmsh.nel);

first_elem = cumsum ([0 hmsh.nel_per_level]) + 1;
last_elem = cumsum ([hmsh.nel_per_level]);
last_dof = cumsum (hspace.ndof_per_level);

if (hspace.ndof == numel(u))
  for ilev = 1:hmsh.nlevels
    if (hmsh.nel_per_level(ilev) > 0)
      msh_level = hmsh.msh_lev{ilev};
      sp_level = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);

      sp_level = change_connectivity_localized_Csub (sp_level, hspace, ilev);

      [errl2_lev, errl2_lev_elem] = ...
        sp_l2_error (sp_level, msh_level, hspace.Csub{ilev}*u(1:last_dof(ilev)), uex);

      errl2 = errl2 + errl2_lev.^2;

      errl2_elem(:,first_elem(ilev):last_elem(ilev))  = errl2_lev_elem;
    end
  end
else
  for ilev = 1:hmsh.nlevels
    if (hmsh.nel_per_level(ilev) > 0)
      msh_level = hmsh.msh_lev{ilev};
      sp_level = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);

      sp_level = change_connectivity_localized_Csub (sp_level, hspace, ilev);
      sp_level = sp_scalar2vector (sp_level, msh_level, 'value', true);

      dofs_to_lev = [];
      for icomp = 1:hmsh.rdim
        dofs_to_lev = union (dofs_to_lev, (icomp-1)*hspace.ndof + (1:last_dof(ilev)));
      end

      u_lev = reshape (u(dofs_to_lev), last_dof(ilev), hmsh.rdim);
      u_lev =  reshape (hspace.Csub{ilev} * u_lev, [], 1);
      [errl2_lev, errl2_lev_elem] = ...
        sp_l2_error (sp_level, msh_level, u_lev, uex);

      errl2 = errl2 + errl2_lev.^2;

      errl2_elem(:,first_elem(ilev):last_elem(ilev))  = errl2_lev_elem;
    end
  end
end

errl2  = sqrt (errl2);

end
