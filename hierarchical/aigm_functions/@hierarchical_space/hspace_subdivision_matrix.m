% HSPACE_SUBDIVISION_MATRIX: compute the matrices for changing basis, from
%                 active functions to B-spline of the tensor product spaces.
%
%   Csub = hspace_subdivision_matrix (hspace, [hmsh], [option])
%
% INPUT:
%
%   hspace:    object representing the hierarchical space (see hierarchical_space)
%   hmsh:      object representing the hierarchical mesh, only needed for the 'reduced' version (see hierarchical_mesh)
%   option:    either 'reduced' (default) or 'full'. The first only uses functions
%               acting on active elements (used for assembly); the second uses the 
%               whole basis of each level (used for plotting), 
%
% OUTPUT:
%
%   Csub:      the matrices for basis change
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
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

function Csub = hspace_subdivision_matrix (hspace, hmsh, option)

if (nargin < 3)
  option = 'reduced';
elseif (~strcmpi (option, 'reduced') && ~strcmpi (option, 'full'))
  warning ('Unknown option. Computing the reduced version')
end

Csub = cell (hspace.nlevels, 1);
Csub{1} = speye (hspace.space_of_level(1).ndof);
Csub{1} = Csub{1}(:,hspace.active{1});

for lev = 2:hspace.nlevels
  I = speye (hspace.space_of_level(lev).ndof);
  aux = matrix_basis_change__ (hspace, lev);
  Csub{lev} = [aux*Csub{lev-1}, I(:,hspace.active{lev})];
  clear aux I
  if (strcmpi (option, 'reduced'))
    functions = sp_get_basis_functions (hspace.space_of_level(lev-1), hmsh.mesh_of_level(lev-1), hmsh.active{lev-1});
    functions = setdiff (1:hspace.space_of_level(lev-1).ndof, functions);
    Csub{lev-1}(functions,:) = 0;
    [i,j,v] = find (Csub{lev-1});
    [m, n] = size (Csub{lev-1});
    Csub{lev-1} = sparse (i, j, v, m, n);
  end
end

if (strcmpi (option, 'reduced'))
  functions = sp_get_basis_functions (hspace.space_of_level(hmsh.nlevels), hmsh.mesh_of_level(hmsh.nlevels), hmsh.active{hmsh.nlevels});
  functions = setdiff (1:hspace.space_of_level(hmsh.nlevels).ndof, functions);
  Csub{hmsh.nlevels}(functions,:) = 0;
  [i,j,v] = find (Csub{hmsh.nlevels});
  [m, n] = size (Csub{hmsh.nlevels});
  Csub{hmsh.nlevels} = sparse (i, j, v, m, n);
end

end
