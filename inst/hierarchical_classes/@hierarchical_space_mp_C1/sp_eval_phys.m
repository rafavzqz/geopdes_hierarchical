% SP_EVAL_PHYS: Compute the value or derivatives of a function from its degrees of freedom, at a given set of points in the physical domain.
%
%   eu = sp_eval_phys (u, hspace, hmsh, geometry, pts, [patch_list], [options]);
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

function eu = sp_eval_phys (u, hspace, hmsh, geometry, pts, patch_list, varargin)

  if (numel (u) ~= hspace.ndof && numel(u) ~= hmsh.rdim*hspace.ndof)
    error ('The number of degrees of freedom of the vector and the space do not match')
  end
  if (numel(geometry) ~= hspace.space_of_level(1).npatch)
    error ('The number of patches of the space and geometry do not coincide')
  end    
  if (~isfield (geometry, 'nurbs'))
    error ('Only implemented for NURBS geometries')
  end

  npts = size (pts, 2);
  ndim = hmsh.ndim;
  rdim = hmsh.rdim;
  npatch = numel (geometry);
  pts_param = zeros (ndim, npts); flag = false (npts, 1);
  
% This is necessary because the NURBS toolbox works in the 3D space
  if (size(pts, 1) < 3)
    pts(end+1:3,:) = 0;
  end

  if (nargin < 6 || isempty (patch_list))
    patch_list = zeros (1, npts);
    for ipnt = 1:npts
      for iptc = 1:npatch
        [pts_param(:,ipnt), flag(ipnt)] = nrbinverse (geometry(iptc).nurbs, pts(:,ipnt));
        if (flag(ipnt))
          patch_list(ipnt) = iptc;
          break
        end
      end
    end
  else
    for ipnt = 1:npts
      iptc = patch_list(ipnt);
      [pts_param(:,ipnt), flag(ipnt)] = nrbinverse (geometry(iptc).nurbs, pts(:,ipnt));
    end
  end

  pts_on_patch = cell (1, npatch);
  msh_aux = struct ('ndim', ndim, 'rdim', rdim, 'boundary', []);
  res = cell (npts, 1);
  for ipnt = 1:npts
    iptc = patch_list(ipnt);
    coords = pts_param(:,ipnt);
    coords_cell = mat2cell (coords, ones(numel(coords), 1));
    for ilev = 1:hmsh.nlevels
      breaks = hmsh.mesh_of_level(ilev).msh_patch{iptc}.breaks;
      for idim = 1:ndim
        ind_elem{idim} = findspan (numel(breaks{idim})-2, 0, coords(idim), breaks{idim}) + 1;
      end
      elem_on_patch = sub2ind (hmsh.mesh_of_level(ilev).msh_patch{iptc}.nel_dir,ind_elem{:});
      elem_number = sum(hmsh.mesh_of_level(ilev).nel_per_patch(1:iptc-1)) + elem_on_patch;

      if (ismember (elem_number, hmsh.active{ilev}))
        Csub = repmat (hspace.Csub(ilev), 1, rdim);
        Csub = blkdiag (Csub{:});
        inds = [];
        ndof_until_lev = sum(hspace.ndof_per_level(1:ilev));
        for idim = 1:rdim
          inds = union (inds, (1:ndof_until_lev)+(idim-1)*hspace.ndof);
        end
        % u_rows = Csub * u(1:sum(hspace.ndof_per_level(1:ilev)*rdim));
        u_rows = Csub * u(inds);

        Csub_rows = [];
        ndof_lev = hspace.space_of_level(ilev).ndof;
        for icomp = 1:rdim
          Csub_rows = union (Csub_rows, (icomp-1)*ndof_lev + hspace.Csub_row_indices{ilev});
        end
        [Cpatch, Cpatch_cols] = sp_compute_Cpatch_vector (hspace.space_of_level(ilev), iptc, rdim);
        [~,iCp,iCs] = intersect (Cpatch_cols, Csub_rows);
        u_ptc_lev = Cpatch(:,iCp) * u_rows(iCs);
          
        sp_vec = sp_vector (repmat (hspace.space_of_level(ilev).sp_patch(iptc), rdim, 1), msh_aux);
        res{ipnt} = sp_eval (u_ptc_lev, sp_vec, geometry(iptc), coords_cell(:)', varargin{:});
        break
      end
    end
  end
  eu = cat(2,res{:});

  % Csub = hspace_subdivision_matrix (hspace, hmsh, 'full');
  % Csub = Csub{end};
  % if (numel(u) == hspace.ndof)
  %   eu = sp_eval_phys (Csub * u, hspace.space_of_level(end), geometry, pts, patch_list, varargin{:});
  % else
  %   Csub = repmat({Csub}, rdim, 1);
  %   eu = sp_eval_phys (blkdiag(Csub{:}) * u, hspace.space_of_level(end), geometry, pts, patch_list, varargin{:});
  % end
end
  
