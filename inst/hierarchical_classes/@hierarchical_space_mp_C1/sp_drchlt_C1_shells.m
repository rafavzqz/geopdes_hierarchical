% SP_DRCHLT_C1_SHELLS: 
%
% USAGE:
%
%  [u_drchlt, drchlt_dofs, kernel_info] = sp_drchlt_C1_shells (hspace, hmsh, refs, components)
%
% INPUT:
%
%  hspace:     object representing the space of hierarchical splines (see hierarchical_space_mp_C1)
%  hmsh:       object representing the hierarchical mesh (see hierarchical_mesh_mp)
%  refs:       reference of the boundaries on which to apply the condition
%  components: vector components on which to apply the zero displacement condition, for each boundary
%
% OUTPUT:
%
%  u:        the computed degrees of freedom
%
% Copyright (C) 2023, Rafael Vazquez
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

function [u_drchlt, drchlt_dofs, kernel_info] = sp_drchlt_C1_shells (hspace, hmsh, refs, drchlt_components)

drchlt_dofs = [];

shifting_index = cumsum ([0 hspace.ndof_per_level]);
boundaries = hmsh.mesh_of_level(1).boundaries;
for iref = 1:numel(refs)
%   href = @(varargin) h(varargin{:}, iref);
  if (~exist('drchlt_components','var'))
    components = 1:3;
  else
    components = drchlt_components{iref};
  end
  bnd_ref = refs(iref);
  scalar_dofs_on_ref = [];
  for bnd_side = 1:boundaries(bnd_ref).nsides
    iptc = boundaries(bnd_ref).patches(bnd_side);
    iside = boundaries(bnd_ref).faces(bnd_side);

    for ilev = 1:hmsh.nlevels
      side_dofs = hspace.space_of_level(ilev).sp_patch{iptc}.boundary(iside).dofs;
      [Cpatch, Cpatch_cols] = sp_compute_Cpatch (hspace.space_of_level(ilev), iptc);

      [~,scalar_dofs] = find (Cpatch(side_dofs,:));
      [~,~,dofs_on_lev] = intersect (Cpatch_cols(scalar_dofs), hspace.active{ilev});
      dofs_on_lev = shifting_index(ilev) + dofs_on_lev;
      scalar_dofs_on_ref = union (scalar_dofs_on_ref, dofs_on_lev);
    end
  end
  for icomp = components
    drchlt_dofs = union (drchlt_dofs, (icomp-1)*hspace.ndof + scalar_dofs_on_ref);
  end
end

dofs_to_remove = [];
vertices_numbers = [];
row_indices = [];
count_vert = 0;
count_fun = 0;

% Check the kernel of vertex functions on Dirichlet boundary vertices
% Pick up the basis function with the max. abs. coeff in the kernel, 
%  remove it from drchlt_dofs, and add the function in the kernel into the
%  internal part (it goes in the output)
B_change_local = [];
n_boundaries = numel(boundaries); % number of boundary edges
global_refs = numel(hspace.space_of_level(1).interfaces) - n_boundaries + refs; % global numbering of Dirichlet boundary edges

vertices = hspace.space_of_level(1).vertices;
for iv = 1:numel(vertices)
  for lev = 1:hspace.nlevels
    if (ismember(hspace.space_of_level(lev).dofs_on_vertex{iv}(1), hspace.active{lev}))
      vertex_level = lev;
      continue
    end
  end
  
  % Loop just over Dirichlet boundary vertices
  if ~isempty(intersect(global_refs, vertices(iv).edges))
    if (vertices(iv).boundary_vertex)
      patches = vertices(iv).patches([1 end]);

      operations = vertices(iv).patch_reorientation([1 end], :);
      indices_loc_R = hspace.space_of_level(vertex_level).sp_patch{patches(1)}.ndof_dir;
      indices_loc_L = hspace.space_of_level(vertex_level).sp_patch{patches(2)}.ndof_dir;
      ndof_dir_p1 = indices_loc_R;
      indices_loc_R = indices_reorientation(indices_loc_R, operations(1, :));
      indices_loc_L = indices_reorientation(indices_loc_L, operations(2, :));
      indices_loc_R = indices_loc_R(:);
      indices_loc_L = indices_loc_L(:);

      Cpatch_ind_R = indices_loc_R([1 2 3]);
      Cpatch_ind_L = indices_loc_L([ndof_dir_p1(1)+1 2*ndof_dir_p1(1)+1]);
%       Cpatch_ind_R = indices_loc_R([1 2 3 space.sp_patch{patches(1)}.ndof_dir(1)+[1 2]]);
%       if (space.vertices(iv).valence_p == 2)
%         Cpatch_ind_L = indices_loc_L([space.sp_patch{patches(2)}.ndof_dir(1)+[1 2] 2*space.sp_patch{patches(2)}.ndof_dir(1)+1]);
%       else
%         Cpatch_ind_L = indices_loc_L([2 space.sp_patch{patches(2)}.ndof_dir(1)+[1 2] 2*space.sp_patch{patches(2)}.ndof_dir(1)+1]);
%       end

      [Cpatch1, Cpatch_cols1] = sp_compute_Cpatch (hspace.space_of_level(vertex_level), patches(1));
      [Cpatch2, Cpatch_cols2] = sp_compute_Cpatch (hspace.space_of_level(vertex_level), patches(2));

      [~,~,inds1] = intersect (hspace.space_of_level(vertex_level).dofs_on_vertex{iv}, Cpatch_cols1);
      [~,~,inds2] = intersect (hspace.space_of_level(vertex_level).dofs_on_vertex{iv}, Cpatch_cols2);

      M_ker = [Cpatch1(Cpatch_ind_R, inds1); ...
               Cpatch2(Cpatch_ind_L, inds2)];

      ker = null(full(M_ker));
      if (~isempty(ker))
        nfun = size(ker,2);
        [~, ind] = max(abs(ker)); % TODO: NOT A GOOD CHOICE (it may be repeated)

        row_inds = count_vert*6 + (1:6);
        B_change_local = blkdiag (B_change_local, ker);

%         dofs_on_vertex = space.dofs_on_vertex{iv};
        dofs_on_vertex = hspace.space_of_level(vertex_level).dofs_on_vertex{iv};
        [~,~,position] = intersect (dofs_on_vertex, hspace.active{vertex_level});
        dofs_on_vertex = shifting_index(vertex_level) + position;

        dofs_to_remove(count_fun+(1:nfun)) = dofs_on_vertex(ind);
        row_indices(row_inds) = dofs_on_vertex;
        vertices_numbers(count_fun+(1:nfun)) = iv;

        count_vert = count_vert + 1;
        count_fun  = count_fun + nfun;
      end
    end
  end
end

kernel_info = struct ('vertices_numbers', vertices_numbers, 'all_vertex_dofs', row_indices, 'quasi_interior_dofs', dofs_to_remove, 'B_change_local', sparse(B_change_local));

dofs_to_remove = [dofs_to_remove(:), dofs_to_remove(:)+hspace.ndof, dofs_to_remove(:)+2*hspace.ndof];
drchlt_dofs = setdiff(drchlt_dofs, dofs_to_remove);

u_drchlt = zeros (numel(drchlt_dofs), 1);

end

function indices = indices_reorientation (ndof_dir, operations)
  ndof = prod (ndof_dir);
  indices = reshape (1:ndof, ndof_dir);
  if (operations(1))
    indices = flipud (indices);
  end
  if (operations(2))
    indices = fliplr (indices);
  end
  if (operations(3))
    indices = indices.';
  end   
end
