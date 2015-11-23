% HSPACE_EVAL_HMSH: Compute the value or the derivatives of a hierarchical spline function, 
%  given by its degrees of freedom, at the points of the corresponding hierarchical mesh.
%
%   [eu, F] = hspace_eval_hmsh (u, hspace, hmsh, [option]);
%   [eu, F] = hspace_eval_hmsh (u, hspace, hmsh, [option]);
%
% INPUT:
%     
%     u:         vector of dof weights
%     hspace:    object defining the discrete space (see hierarchical_space)
%     hmsh:      object representing the hierarchical mesh (see hierarchical_mesh)
%     option:    accepted options are 'value' (default), 'gradient', 'laplacian'
%
% OUTPUT:
%
%     eu: the function evaluated at the points in the hierarchical mesh
%     F:  grid points in the physical domain, that is, the mapped points
% 
% Copyright (C) 2015 Rafael Vazquez
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

function [eu, F] = hspace_eval_hmsh (u, hspace, hmsh, varargin)

  if (hspace.ncomp ~= 1)
    error ('hspace_eval_hmsh: Not implemented for vector valued spaces')
  end

  if (nargin == 3)
    option = 'value';
  else
    option = varargin{1};
  end

  value = false; gradient = false; laplacian = false;
  switch (lower (option))
    case {'value'}
      eval_fun = @(U, SP, MSH) sp_eval_msh (U, SP, MSH, 'value');
      catdir = 2;
      value = true;
    case {'gradient'}
      eval_fun = @(U, SP, MSH) sp_eval_msh (U, SP, MSH, 'gradient');
      catdir = 3;
      gradient = true;
    case {'laplacian'}
      eval_fun = @(U, SP, MSH) sp_eval_msh (U, SP, MSH, 'laplacian');
      catdir = 2;
      laplacian = true;
%     case {'divergence'}
%     case {'curl'}
    otherwise
      error ('hspace_eval_hmsh: unknown option to evaluate')
  end

  F = []; eu = [];

  last_dof = cumsum (hspace.ndof_per_level);
  for ilev = 1:hmsh.nlevels % Active levels
    if (hmsh.nel_per_level(ilev) > 0)
      msh_level = hmsh.msh_lev{ilev};
%      sp_level = sp_evaluate_element_list (hspace.space_of_level(ilev), msh_level, 'value', value, 'gradient', gradient);
      sp_level = sp_evaluate_element_list (hspace.space_of_level(ilev), msh_level, 'value', value, 'gradient', gradient, 'laplacian', laplacian);
        
      [eu_lev, F_lev] = eval_fun (hspace.C{ilev}*u(1:last_dof(ilev)), sp_level, msh_level);
      eu = cat (catdir, eu, eu_lev);
      F = cat (3, F, F_lev);
    end
  end
  
end