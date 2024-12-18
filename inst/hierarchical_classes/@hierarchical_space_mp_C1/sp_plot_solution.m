% SP_PLOT_SOLUTION: Plot the computed solution, given the degrees of freedom.
%
%   [eu, F] = sp_plot_solution (u, space, geometry, hmsh, pts);
%   [eu, F] = sp_plot_solution (u, space, geometry, hmsh, [npts]);
%
% INPUT:
%     
%     u:           vector of dof weights
%     hspace:      object defining the discrete space (see hierarchical_space_mp_C1)
%     geometry:    geometry structure (see mp_geo_load)
%     hmsh:        object defining the hierarchical mesh (see hierarchical_mesh_mp)
%     pts:         cell array with coordinates of points along each parametric direction
%     npts:        number of points along each parametric direction
%
%    This function only plots the value of the solution. To plot other
%     quantities, such as the gradient, compute them with sp_eval.
%
%    The current version is very unefficient, as it passes the solution to
%     the tensor-product space of the finest level.
%
% Copyright (C) 2015-2024 Rafael Vazquez
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

function sp_plot_solution (u, hspace, geometry, hmsh, varargin)

if (numel(u) == hspace.ndof)
  sp_lev = hspace.space_of_level(hspace.nlevels);
  C = hspace_subdivision_matrix (hspace, hmsh, 'full'); % FIX: remove hmsh
  u_lev =  C{hspace.nlevels} * u;
  sp_plot_solution (u_lev, sp_lev, geometry, varargin{:});
else
  error ('Not implemented yet')
end

end