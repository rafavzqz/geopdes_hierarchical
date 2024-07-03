% ADAPTIVITY_ESTIMATE_LAPLACE_H_H2: compute the estimator by solving the problem in a globally refined mesh.
%  The implementation is unefficient, but the estimator should work properly.
%
% USAGE:
%
%   est = adaptivity_estimate_laplace_h_h2 (u, hmsh, hspace, problem_data, method_data)
%
% INPUT:
%
%   u:            degrees of freedom
%   hmsh:         object representing the hierarchical mesh (see hierarchical_mesh_mp)
%   hspace:       object representing the space of hierarchical splines (see hierarchical_space_mp_C1)
%   problem_data: a structure with data of the problem (see adaptivity_solve_laplace_mp_C1).
%   method_data:  a structure with discretization data (only used for penalization in Nitsche's method)
%
% OUTPUT:
%
%   est: computed a posteriori error indicator
%
%
% Copyright (C) 2023 Rafael Vazquez
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


function est = adaptivity_estimate_laplace_h_h2 (u, hmsh, hspace, problem_data, method_data)

if (~isempty (hmsh.active{end}))
  hmsh = hmsh_add_new_level (hmsh);
  hspace = hspace_add_new_level (hspace, hmsh);
end

adaptivity_data.flag = 'elements';
[hmsh_h2, hspace_h2, Cref] = adaptivity_refine (hmsh, hspace, hmsh.active, adaptivity_data);

u_h2 = adaptivity_solve_laplace_mp_C1 (hmsh_h2, hspace_h2, problem_data, method_data);

% Compute the estimator on the fine mesh, and then pass to the coarse mesh
zeroex = @(varargin) zeros (size(varargin{1}));
zerogex = @(varargin) zeros ([hmsh.rdim, size(varargin{1})]);
[~,~,~,~,~,est_elems_h2] = sp_h1_error (hspace_h2, hmsh_h2, u_h2 - Cref*u, zeroex, zerogex);

first_elem = cumsum ([0 hmsh.nel_per_level]) + 1;
last_elem = cumsum ([hmsh.nel_per_level]);
last_elem_h2 = cumsum ([hmsh_h2.nel_per_level]);
est = zeros (hmsh.nel, 1);
nchildren = prod (hmsh.nsub);
for ilev = 1:hmsh.nlevels-1
  [~,~,children_inds] = hmsh_get_children (hmsh, ilev, hmsh.active{ilev});
  [~,children_pos] = ismember (children_inds, hmsh_h2.active{ilev+1});
  inds_hmsh = first_elem(ilev):last_elem(ilev);
  inds_hmsh_h2 = last_elem_h2(ilev) + children_pos;
  est_h2 = reshape (est_elems_h2(inds_hmsh_h2), nchildren, numel(inds_hmsh));
  est(inds_hmsh) = sqrt (sum (est_h2.^2, 1));
end

end