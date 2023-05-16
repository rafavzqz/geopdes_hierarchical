% SP_EVAL_PHYS: Compute the value or derivatives of a function from its degrees of freedom, at a given set of points in the physical domain.
%
%   eu = sp_eval_phys (u, space, geometry, pts, [options]);
%   [eu, pts_list] = sp_eval_phys (u, space, geometry, pts, [options]);
%
% INPUT:
%     
%     u:           vector of dof weights
%     space:       object defining the discrete space (see sp_scalar)
%     geometry:    geometry structure (see geo_load)
%     pts:         array (rdim x npts) with coordinates of points
%     npts:        number of points along each parametric direction
%     options:     cell array with the fields to plot
%                   accepted options are 'value' (default), 'gradient', 'laplacian', 'hessian'
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

function [eu, pts_list] = sp_eval_phys (u, hspace, hmsh, geometry, pts, options)

  if (numel (u) ~= hspace.ndof)
    error ('The number of degrees of freedom of the vector and the space do not match')
  end
  if (~isfield (geometry, 'nurbs'))
    error ('Only implemented for NURBS geometries')
  end

  if (nargin < 6)
    options = {'value'};
  end
  if (~iscell (options))
    options = {options};
  end

  if (isa (hspace.space_of_level(1), 'sp_scalar'))
    is_scalar = true;
  else
    is_scalar = false;
  end
  
  nopts = numel (options);
  ndim = hmsh.ndim;
  rdim = hmsh.rdim;
  
  nurbs = geometry.nurbs;

  npts = size (pts, 2);
  pts_param = zeros (ndim, npts); flag = false (npts, 1);
% This is necessary because the NURBS toolbox works in the 3D space
  if (size(pts, 1) < 3)
    pts(end+1:3,:) = 0;
  end
  
  for ipt = 1:npts
    [pts_param(:,ipt), flag(ipt)] = nrbinverse (nurbs, pts(:,ipt));
  end
  
  if (~all(flag))
    warning ('Some of the points are not on the geometry')
    pts_param = pts_param(:,flag);
    npts_on_F = size (pts_param, 2);
    if (nargout == 2)
      npts_out = npts_on_F;
      pts_list = find (flag(:)');
      inds_on_F = 1:npts_out;
    else
      npts_out = npts;
      inds_on_F = find (flag);
    end
  else
    npts_on_F = npts;
    npts_out = npts;
    inds_on_F = 1:npts;
    pts_list = 1:npts;
  end
  
  for ilev = 1:hspace.nlevels
    for idim = 1:ndim
      zeta{ilev}{idim} = unique (hmsh.mesh_of_level(ilev).breaks{idim});
      ind{ilev}{idim} = findspan (numel(zeta{ilev}{idim})-2, 0, pts_param(idim,:), zeta{ilev}{idim}) + 1;
    end
  end

  [eu, eunum, eusize, fields_to_plot] = set_output_sizes (ndim, rdim, npts_out, hspace.ncomp, is_scalar, options);
  
  last_dof = cumsum (hspace.ndof_per_level);
  for ilev = 1:hmsh.nlevels
    u_lev = hspace.Csub{ilev}*u(1:last_dof(ilev));
    elem_list = sub2ind (hmsh.mesh_of_level(ilev).nel_dir, ind{ilev}{:});
    [~,pts_on_lev,~] = intersect (elem_list, hmsh.active{ilev});
    for ipt = 1:numel(pts_on_lev)
      point = pts_on_lev(ipt);
      for idim = 1:ndim
        brk{idim} = zeta{ilev}{idim}(ind{ilev}{idim}(point):ind{ilev}{idim}(point)+1);
      end
      msh = msh_cartesian (brk, num2cell(pts_param(:,point)), [], geometry, 'boundary', false);
      sp  = hspace.space_of_level(ilev).constructor (msh);
      
      msh_col = msh_evaluate_element_list (msh, 1);
      if (is_scalar)
        sp_col  = sp_evaluate_element_list (sp, msh_col, 'value', fields_to_plot.value, 'gradient', fields_to_plot.gradient, ...
                'laplacian', fields_to_plot.laplacian, 'hessian', fields_to_plot.hessian);
      else
        sp_col  = sp_evaluate_element_list (sp, msh_col, 'value', fields_to_plot.value, 'gradient', fields_to_plot.gradient, ...
                'curl', fields_to_plot.curl, 'divergence', fields_to_plot.divergence);
      end
      sp_col = change_connectivity_localized_Csub (sp_col, hspace, ilev);
      eu_aux = sp_eval_msh (u_lev, sp_col, msh_col, options);

      for iopt = 1:nopts
        eu{iopt}(eunum{iopt}{:},inds_on_F(point)) = eu_aux{iopt};
      end  
    end
  end
  
  for iopt = 1:nopts
    eu{iopt} = reshape (eu{iopt}, [eusize{iopt}, 1]); % The extra 1 makes things work also in 1D
  end

  if (nopts == 1)
    eu = eu{1};
  end

end

function [eu, eunum, eusize, fields_to_plot] = set_output_sizes (ndim, rdim, npts, ncomp, is_scalar, options)

  nopts = numel(options);
  eu = cell (nopts, 1); eunum = eu; eusize = eu;

  if (is_scalar)
    value = false; grad = false; laplacian = false; hessian = false;
    for iopt = 1:nopts
      switch (lower (options{iopt}))
        case 'value'
          eu{iopt} = NaN (1, npts);
          eunum{iopt} = {1};
          eusize{iopt} = npts;
          value = true;
        case 'gradient'
          eu{iopt} = NaN (rdim, npts);
          eunum{iopt} = {1:rdim};
          eusize{iopt} = [rdim, npts];
          grad = true;
        case 'laplacian'
          eu{iopt} = NaN (1, npts);
          eunum{iopt} = {1};
          eusize{iopt} = npts;
          laplacian = true;
        case 'hessian'
          eu{iopt} = NaN (rdim, rdim, npts);
          eunum{iopt} = {1:rdim, 1:rdim};
          eusize{iopt} = [rdim, rdim, npts];
          hessian = true;
      end
    end
    fields_to_plot = struct ('value', value, 'gradient', grad, 'laplacian', laplacian, 'hessian', hessian);
  else
    value = false; grad = false; curl = false; divergence = false;
    for iopt = 1:numel(options)
      switch (lower (options{iopt}))
        case 'value'
          eu{iopt} = NaN (ncomp, npts);
          eunum{iopt} = {1:ncomp};
          eusize{iopt} = [ncomp, npts];
          value = true;
        case 'gradient'
          eu{iopt} = NaN (ncomp, rdim, npts);
          eunum{iopt} = {1:ncomp, 1:rdim};
          eusize{iopt} = [ncomp, rdim, npts];
          grad = true;
        case 'curl'
          if (ndim == 2 && rdim == 2)
            eu{iopt} = NaN (1, npts);
            eunum{iopt} = {1};
            eusize{iopt} = npts;
          elseif (ndim == 3 && rdim == 3)
            eu{iopt} = NaN (ncomp, rdim, npts);
            eunum{iopt} = {1:ncomp};
            eusize{iopt} = [rdim, npts];
          end
          curl = true;
        case 'divergence'
          eu{iopt} = NaN (1, npts);
          eunum{iopt} = {1};
          eusize{iopt} = npts;
          divergence = true;
      end
    end
    fields_to_plot = struct ('value', value, 'gradient', grad, 'curl', curl, 'divergence', divergence);
  end
end
