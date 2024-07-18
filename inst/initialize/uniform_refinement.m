% UNIFORM_REFINEMENT: construct a hierarchical mesh and space, in the
%   parametric domain, from a given set of either active or refined cells.
%
%  [hmsh, hspace] = uniform_refinement (hmsh, hspace, n_refinements, adaptivity_data)
%
% INPUT:
%
%   hmsh:          hierarchical mesh of one single level
%   hspace:        hierarchical space of one single level
%   n_refinements: times to refine uniformly (number of levels minus one)
%   adaptivity_data: for this function it must contain
%     - max_level: maximum number of allowed levels
%
% OUTPUT:
%
%   hmsh:   hierarchical mesh with all active elements of the finest level
%   hspace: hierarchical space with all active functions of the finest level
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
%

function [hmsh, hspace] = uniform_refinement(hmsh, hspace, n_refinements, adaptivity_data)

% Since the refinement is global, we simplify by marking elements
  adaptivity_data.flag = 'elements';
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
