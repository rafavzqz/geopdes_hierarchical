% SP_EVAL: Compute the value or the derivatives of a hierarchical spline function, given by its degrees of freedom, at a given set of points.
%
%   [eu, F] = sp_eval (u, hspace, geometry, pts, [option]);
%   [eu, F] = sp_eval (u, hspace, geometry, npts, [option]);
%
% INPUT:
%     
%     u:         vector of dof weights
%     hspace:    object defining the discrete space (see hierarchical_space)
%     geometry:  geometry structure (see geo_load)
%     pts:       cell array with coordinates of points along each parametric direction
%     npts:      number of points along each parametric direction
%     option:    accepted options are 'value' (default), 'gradient', 'laplacian'
%
% OUTPUT:
%
%     eu: the function evaluated at the given points 
%     F:  grid points in the physical domain, that is, the mapped points
% 
%    The current version is very unefficient, as it passes the solution to
%     the tensor-product space of the finest level.
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

function [eu, F] = sp_eval (u, hspace, hmsh, geometry, npts, varargin)

  sp_lev = hspace.space_of_level(hspace.nlevels);
  C = hspace_subdivision_matrix (hspace);
  u_lev =  C{hspace.nlevels} * u;
  [eu, F] = sp_eval (u_lev, sp_lev, geometry, npts, varargin{:});

% output_cell = true;  
% if (nargin == 5)
%     varargin = {'value'};
%     output_cell = false;  
% elseif (~iscell (varargin))
%     varargin = {varargin};
%     output_cell = false;
% end
% 
% nopts = numel (varargin);
% 
% % For vector-valued spaces, the value of catdir is then corrected by adding one
% value = false; gradient = false; laplacian = false;
% hessian = false; curl = false; divergence = false;
% for iopt = 1:nopts
%     switch (lower (varargin{iopt}))
%       case 'value'
%         value = true;
%         catdir(iopt) = 2;
%       case 'gradient'
%         gradient = true;
%         catdir(iopt) = 3;
%       case 'laplacian' % Only for scalars, at least for now
%         laplacian = true;
%         catdir(iopt) = 2;
%       case 'hessian'
%         hessian = true;
%         catdir(iopt) = 4;
%       case 'curl' % Only for vectors
%         curl = true;
%         if (hspace.space_of_level(1).ncomp_param == 2)
%           catdir(iopt) = 1;
%         elseif (hspace.space_of_level(1).ncomp_param == 3)
%           catdir(iopt) = 2;
%         end
%       case 'divergence' % Only for vectors
%         divergence = true;
%         catdir(iopt) = 1;
%       otherwise
%         error ('hspace_eval_hmsh: unknown option: %s', varargin{iopt})
%     end
% end
% if (hspace.ncomp ~= 1)
%     catdir = catdir + 1;
% end
% 
% F = []; eu = cell (nopts, 1);
% 
% last_dof = cumsum (hspace.ndof_per_level);
% 
% for iLev = 1:hspace.nlevels
%   if (hspace.ndof_per_level(iLev) > 0)  
%     sp_level = hspace.space_of_level(iLev);
% 
% %     indices = unique (sp_level.connectivity);
% %     [~,position] = ismember (sp_level.connectivity, indices);
%     fun_on_active = sp_get_basis_functions (hspace.space_of_level(iLev), hmsh.mesh_of_level(iLev), hmsh.active{iLev});
%     fun_on_deact = sp_get_basis_functions (hspace.space_of_level(iLev), hmsh.mesh_of_level(iLev), hmsh.deactivated{iLev});
%     fun_on_deact = union (fun_on_active, fun_on_deact);
%     sp_level.ndof = numel (fun_on_deact);
% %     sp_level.connectivity = position;
%     
%     [eu_lev, F_lev] = sp_eval (hspace.Csub{iLev} * u(1:last_dof(iLev)), sp_level, geometry, npts, varargin{:});
%     if (nopts == 1)
%         eu_lev = cell (nopts, 1);
%         eu_lev{1} = eu_lev;
%     end
%   
%     for iopt = 1:nopts
%         eu{iopt} = cat (catdir(iopt), eu{iopt}, eu_lev{iopt});
%     end
%     F = cat (3, F, F_lev);
%   end
% end

end