% SP_EVAL_PHYS: Compute the value or derivatives of a function from its degrees of freedom, at a given set of points in the physical domain.
%
%   eu = sp_eval_phys (u, space, geometry, pts, [patch_list], [options]);
%
% INPUT:
%     
%     u:           vector of dof weights
%     hspace:      object defining the discrete space (see hierarchical_space_mp)
%     hmsh:        object defining the hierarchical mesh (see hierarchical_mesh_mp)
%     geometry:    geometry structure (see mp_geo_load)
%     pts:         array (rdim x npts) with coordinates of points
%     patch_list:  patch to which each point begins. If empty, the patch
%                   will be found using nrbinverse
%     options:     cell array with the fields to plot
%                   accepted options for scalars are 'value' (default), 'gradient', 'laplacian', 'hessian'
%                   accepted options for vectors are 'value' (default), 'gradient', 'curl', 'divergence'
%
% OUTPUT:
%
%     eu:       cell-array with the fields evaluated at the given points.
%     pts_list: list of points that are on the geometry, for which the
%                values are computed.
%  If there is only one output argument, points not on the geometry get a NaN value.
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

function eu = sp_eval_phys (u, hspace, hmsh, geometry, pts, varargin)

  if (numel (u) ~= hspace.ndof && numel(u) ~= hmsh.rdim*hspace.ndof)
    error ('The number of degrees of freedom of the vector and the space do not match')
  end
  if (numel(geometry) ~= hspace.space_of_level(1).npatch)
    error ('The number of patches of the space and geometry do not coincide')
  end    
  if (~isfield (geometry, 'nurbs'))
    error ('Only implemented for NURBS geometries')
  end

  Csub = hspace_subdivision_matrix (hspace, hmsh, 'full');
  Csub = Csub{end};
  if (numel(u) == hspace.ndof)
    eu = sp_eval_phys (Csub * u, hspace.space_of_level(end), geometry, pts, varargin{:});
  else
    Csub = repmat({Csub}, hmsh.rdim, 1);
    eu = sp_eval_phys (blkdiag(Csub{:}) * u, hspace.space_of_level(end), geometry, pts, varargin{:});
  end
  
end
  
