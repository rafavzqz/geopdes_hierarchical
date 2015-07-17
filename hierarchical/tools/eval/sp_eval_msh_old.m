% SP_EVAL_MSH: Evaluate a function, given by its degrees of freedom, at the points given by a msh object.
%
%   [eu, F] = sp_eval_msh (u, space, msh);
%
% INPUT:
%     
%     u:         vector of dof weights
%     space:     object defining the discrete space (see sp_bspline)
%     msh:       object defining the points where to evaluate (see msh_cartesian)
%
% OUTPUT:
%
%     eu: the function evaluated in the given points 
%     F:  grid points in the physical domain, that is, the mapped points
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2013, 2015 Rafael Vazquez
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

function [eu, F] = sp_eval_msh_old (u, space, msh)

  F  = zeros (msh.rdim, msh.nqn, msh.nel);
  eu = zeros (space.ncomp, msh.nqn, msh.nel);

  uc = zeros (size (space.connectivity));
  uc(space.connectivity~=0) = ...
      u(space.connectivity(space.connectivity~=0));
  weight = reshape (uc, [1, 1, space.nsh_max, msh.nel]);

  space.shape_functions = reshape (space.shape_functions, space.ncomp, ...
                                      msh.nqn, space.nsh_max, msh.nel);

  F(:,:,:) = msh.geo_map;
  eu(:,:,:) = reshape (sum (bsxfun (@times, weight, ...
         space.shape_functions), 3), space.ncomp, msh.nqn, msh.nel);

  if (space.ncomp == 1)
    eu = reshape (eu, msh.nqn, msh.nel);
  end

end
