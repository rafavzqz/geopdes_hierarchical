% SP_H2_EQUIV_LAP_ERROR: Evaluate the error in H^2 equivalent seminorm with laplacian, H^1 and L^2 norms.
%
%   [errh2, errh1, errl2, errh2s, errh1s, errh2_elem, errh1_elem, errl2_elem, errh2s_elem, errh1s_elem] = ...
%     sp_h2_equiv_lap_error (hspace, hmsh, u, uex, graduex, lapuex)
%
% INPUT:
%
%    hspace:  object defining the discrete space (see hierarchical_space)
%    hmsh:    object representing the hierarchical mesh (see hierarchical_mesh)
%    u:       vector of dof weights
%    uex:     function handle to evaluate the exact solution
%    graduex: function handle to evaluate the gradient of the exact solution
%    lapuex:  function handle to evaluate the laplacian of the exact solution
%
% OUTPUT:
%
%     errh2:       error in (equivalent) H^2 norm
%     errh1:       error in H^1 norm
%     errl2:       error in L^2 norm
%     errh2s:      error in H^2 seminorm
%     errh1s:      error in H^1 seminorm
%     errh2_elem:  error in (equivalent) H^2 norm, for each single element
%     errh1_elem:  error in H^1 norm, for each single element
%     errl2_elem:  error in L^2 norm, for each single element
%     errh2s_elem: error in H^2 seminorm, for each single element
%     errh1s_elem: error in H^1 seminorm, for each single element
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2023 Rafael Vazquez
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

function [errh2, errh1, errl2, errh2s, errh1s, errh2_elem, errh1_elem, errl2_elem, errh2s_elem, errh1s_elem] = sp_h2_equiv_lap_error (hspace, hmsh, u, uex, graduex, lapuex)

if (numel(u) ~= hspace.ndof)
  error ('Wrong size of the vector of degrees of freedom')
end

errh2s = 0; errh2 = 0; errh1 = 0; errl2 = 0; errh1s = 0;
errh2_elem = zeros (1, hmsh.nel); errh1_elem = zeros (1, hmsh.nel); errl2_elem = zeros (1, hmsh.nel); errh2s_elem = zeros (1, hmsh.nel); errh1s_elem = zeros (1, hmsh.nel);

first_elem = cumsum ([0 hmsh.nel_per_level]) + 1;
last_elem = cumsum ([hmsh.nel_per_level]);
last_dof = cumsum (hspace.ndof_per_level);
for ilev = 1:hmsh.nlevels
  if (hmsh.nel_per_level(ilev) > 0)
    msh_level = hmsh.msh_lev{ilev};
    sp_level = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true, 'gradient', true, 'laplacian', true);

    sp_level = change_connectivity_localized_Csub (sp_level, hspace, ilev);

    [errh2_lev, errh1_lev, errl2_lev, errh2s_lev, errh1s_lev, errh2_lev_elem, errh1_lev_elem, errl2_lev_elem, errh2s_lev_elem, errh1s_lev_elem] = ...
      sp_h2_equiv_lap_error (sp_level, msh_level, hspace.Csub{ilev}*u(1:last_dof(ilev)), uex, graduex, lapuex);

    errh2 = errh2 + errh2_lev.^2;
    errh1 = errh1 + errh1_lev.^2;
    errl2 = errl2 + errl2_lev.^2;
    errh1s = errh1s + errh1s_lev.^2;
    errh2s = errh2s + errh2s_lev.^2;

    errh2_elem(:,first_elem(ilev):last_elem(ilev))  = errh2_lev_elem;    
    errh1_elem(:,first_elem(ilev):last_elem(ilev))  = errh1_lev_elem;
    errl2_elem(:,first_elem(ilev):last_elem(ilev))  = errl2_lev_elem;
    errh1s_elem(:,first_elem(ilev):last_elem(ilev)) = errh1s_lev_elem;
    errh2s_elem(:,first_elem(ilev):last_elem(ilev)) = errh2s_lev_elem;

  end
end
errh2  = sqrt (errh2);
errh2s  = sqrt (errh2s);
errh1  = sqrt (errh1);
errl2  = sqrt (errl2);
errh1s = sqrt (errh1s);

end
