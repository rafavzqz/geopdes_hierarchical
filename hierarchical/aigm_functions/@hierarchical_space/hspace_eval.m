% HSPACE_EVAL: Compute the value or the derivatives of a hierarchical spline function, given by its degrees of freedom, at a given set of points.
%
%   [eu, F] = hspace_eval (u, hspace, geometry, pts, [option]);
%   [eu, F] = hspace_eval (u, hspace, geometry, npts, [option]);
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

function [eu, F] = hspace_eval (u, hspace, geometry, npts, varargin)

  if (hspace.ncomp ~= 1)
    error ('hspace_eval: Not implemented for vector valued spaces')
  end

  if (nargin == 4)
    option = 'value';
  else
    option = varargin{1};
  end

  ndim = numel (npts);
  
% Temporary solution, to be fixed using "isprop" after defining the
%  classes with classdef
  sp_coarse = hspace.space_of_level(1);
  
  knt = cell (ndim, 1);
  if (isfield (struct(sp_coarse), 'knots'))
    for idim=1:ndim
      knt{idim} = sp_coarse.knots{idim}(sp_coarse.degree(idim)+1:end-sp_coarse.degree(idim));
    end
%   elseif (isfield (struct(space), 'scalar_spaces'))
%     for idim=1:ndim
%       knt{idim} = sp_coarse.scalar_spaces{1}.knots{idim}(sp_coarse.scalar_spaces{1}.degree(idim)+1:end-sp_coarse.scalar_spaces{1}.degree(idim));
%     end
  else
    for idim=1:ndim; knt{idim} = [0 1]; end
  end
  
  if (iscell (npts))
    pts = npts;
    npts = cellfun (@numel, pts);
  elseif (isvector (npts))
    for idim = 1:ndim
      pts{idim} = linspace (knt{idim}(1), knt{idim}(end), npts(idim));
    end
  end

  for jj = 1:ndim
    pts{jj} = pts{jj}(:)';
    if (numel (pts{jj}) > 1)
      brk{jj} = [knt{jj}(1), pts{jj}(1:end-1) + diff(pts{jj})/2, knt{jj}(end)];
    else
      brk{jj} = [knt{jj}(1) knt{jj}(end)];
    end
  end

  msh = msh_cartesian (brk, pts, [], geometry, 'boundary', false);
  
  first_dof = cumsum ([0 hspace.ndof_per_level]) + 1;

  switch (lower (option))
    case {'value'}
      eu = squeeze (zeros ([hspace.ncomp, npts]));
      eval_fun = @(U, SP) sp_eval_msh (U, SP, msh);
    case {'gradient'}
      eu = squeeze (zeros ([hspace.ncomp, msh.rdim, npts]));
      eval_fun = @(U, SP) sp_eval_grad_msh (U, SP, msh);
    case {'laplacian'}
      eu = squeeze (zeros ([hspace.ncomp, npts]));
      eval_fun = @(U, SP) sp_eval_lapl_msh (U, SP, msh);
%     case {'divergence'}
%     case {'curl'}
    otherwise
      error ('hspace_eval: unknown option to evaluate')
  end

  for ilev = 1:hspace.nlevels
    sp_lev = hspace.space_of_level(ilev).constructor (msh);
    u_lev  = zeros (sp_lev.ndof, 1);
    u_lev(hspace.active{ilev}) = u(first_dof(ilev):first_dof(ilev+1)-1);
    
    [eu_lev, F] = eval_fun (u_lev, sp_lev);
    eu = eu + reshape (eu_lev, size(eu));
  end
  F  = reshape (F, [msh.rdim, npts]);
