% HSPACE_GET_CHILDREN: compute the children of a given set of functions of the same level.
%
%     [children, flag, children_of_function] = hspace_get_children (hspace, lev, ind)
%
% Get the children of a given set of basis functions of level lev, with the
%  subdivision given by the "Proj" matrices
% All the children functions are stored in the same array.
%
% INPUT:
%
%     hspace: the hierarchical space (see hierarchical_space_mp_C1)
%     lev:    level of the functions to refine
%     ind:    indices of the functions in the multipatch space of level lev
%
% OUTPUT:
%
%     children: indices of the children, with the numbering of the tensor product space
%     flag:     a flag to tell whether all the input functions are active (1) 
%               active or deactivated (2), or if there is any passive function (0)
%     children_of_function: cell array with the children of each function
%
% Copyright (C) 2015, 2016, 2017 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2018--2022 Cesare Bracco, Rafael Vazquez
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

function [children, flag, children_of_function] = hspace_get_children (hspace, lev, ind)

ndim = size (hspace.Proj{1}, 2);
npatch = hspace.space_of_level(1).npatch;
interfaces = hspace.space_of_level(1).interfaces;
nint = numel (interfaces);
vertices = hspace.space_of_level(1).vertices;
nvert = numel (vertices);

children = [];
children_of_function = cell (numel(ind), 1);
ind_sub = cell (ndim, 1);

% Children of interior functions
aux = cell (ndim, 1);
for iptc = 1:npatch
  [~,int_patch_dofs_c] = sp_get_functions_on_patch (hspace.space_of_level(lev), iptc);
%   int_patch_dofs_c = hspace.space_of_level(lev).dofs_on_patch{iptc};
  [~,indices,position] = intersect (int_patch_dofs_c, ind);
  
  indices_tp = sp_get_local_interior_functions (hspace.space_of_level(lev), iptc);
  indices_tp = indices_tp(indices);
%   indices_tp = hspace.space_of_level(lev).interior_dofs_per_patch{iptc}(indices);

  [ind_sub{:}] = ind2sub ([hspace.space_of_level(lev).sp_patch{iptc}.ndof_dir, 1], indices_tp); % The extra 1 makes it work in any dimension

  for ii = 1:numel(ind_sub{1})
    for idim = 1:ndim
      aux{idim} = find (hspace.Proj{lev, iptc}{idim}(:,ind_sub{idim}(ii)));
    end
    [z{1:ndim}] = ndgrid (aux{:});
    auxI = sub2ind ([hspace.space_of_level(lev+1).sp_patch{iptc}.ndof_dir, 1], z{:});
    local_interior_dofs = sp_get_local_interior_functions (hspace.space_of_level(lev+1), iptc);
    [~,local_indices,~] = intersect (local_interior_dofs, auxI);
%     [~,local_indices,~] = intersect (hspace.space_of_level(lev+1).interior_dofs_per_patch{iptc}, auxI);
    [~,interior_dofs_patch] = sp_get_functions_on_patch (hspace.space_of_level(lev+1), iptc);
    children_of_this_function = interior_dofs_patch(local_indices(:));
%     children_of_this_function = hspace.space_of_level(lev+1).dofs_on_patch{iptc}(local_indices(:));
    children = union (children, children_of_this_function);
    children_of_function{position(ii)} = children_of_this_function(:).';
  end
end

% Children of edge functions
for iint = 1:nint
  patches = [interfaces(iint).patch1 interfaces(iint).patch2];
  sides = [interfaces(iint).side1 interfaces(iint).side2];
  
  [ind_on_edge,~,position] = intersect (hspace.space_of_level(lev).dofs_on_edge{iint}, ind);
  
  possible_children = [];
  ind_sub = cell (ndim, 1);
  for iptc = 1:numel(patches)
    Proj_patch = hspace.Proj{lev, patches(iptc)};
    ndof_dir_bsp_c = hspace.space_of_level(lev).sp_patch{patches(iptc)}.ndof_dir;
    ndof_dir_bsp_f = hspace.space_of_level(lev+1).sp_patch{patches(iptc)}.ndof_dir;
    if (sides(iptc) == 1)
      ind_sub{1} = [1 2];
      ind_sub{2} = 1:ndof_dir_bsp_c(2);
    elseif (sides(iptc) == 2)
      ind_sub{1} = [ndof_dir_bsp_c(1)-1 ndof_dir_bsp_c(1)];
      ind_sub{2} = 1:ndof_dir_bsp_c(2);
    elseif (sides(iptc) == 3)
      ind_sub{1} = 1:ndof_dir_bsp_c(1);
      ind_sub{2} = [1 2];
    elseif (sides(iptc) == 4)
      ind_sub{1} = 1:ndof_dir_bsp_c(1);
      ind_sub{2} = [ndof_dir_bsp_c(2)-1 ndof_dir_bsp_c(2)];
    end
    for idim = 1:ndim
      [aux{idim},~] = find (Proj_patch{idim}(:,ind_sub{idim}));
      aux{idim} = unique (aux{idim});
    end
    [z{1:ndim}] = ndgrid (aux{:});
    auxI = sub2ind ([ndof_dir_bsp_f, 1], z{:});
    local_interior_dofs = sp_get_local_interior_functions (hspace.space_of_level(lev+1), patches(iptc));
    [~,local_indices,~] = intersect (local_interior_dofs, auxI);
%     [~,local_indices,~] = intersect (hspace.space_of_level(lev+1).interior_dofs_per_patch{patches(iptc)}, auxI);
    [~,interior_dofs_patch] = sp_get_functions_on_patch (hspace.space_of_level(lev+1), patches(iptc));
    possible_children = union (possible_children, interior_dofs_patch(local_indices(:)));
%     possible_children = union (possible_children, hspace.space_of_level(lev+1).dofs_on_patch{patches(iptc)}(local_indices(:)));
  end

% %   possible_children = cell2mat (hspace.space_of_level(lev+1).dofs_on_patch(patches));
  possible_children = union (possible_children, hspace.space_of_level(lev+1).dofs_on_edge{iint});
  ref_matrix = matrix_basis_change__ (hspace, lev+1, ind_on_edge, possible_children);

  for ii = 1:numel(ind_on_edge)
    [auxI,~]= find(ref_matrix(:,ii));
    children_of_this_function = possible_children(auxI(:));
    children = union (children, children_of_this_function);
    children_of_function{position(ii)} = children_of_this_function(:).';
  end
end

% Children of vertex functions (I use some brute force for interior functions)
for iv = 1:nvert
  patches = vertices(iv).patches;
  edges = vertices(iv).edges;

  [ind_on_vertex,~,position] = intersect (hspace.space_of_level(lev).dofs_on_vertex{iv}, ind);
  
  possible_children = [];
  ind_sub = cell (ndim, 1);
  for iptc = 1:numel(patches)
    Proj_patch = hspace.Proj{lev, patches(iptc)};
    ndof_dir_bsp_c = hspace.space_of_level(lev).sp_patch{patches(iptc)}.ndof_dir;
    ndof_dir_bsp_f = hspace.space_of_level(lev+1).sp_patch{patches(iptc)}.ndof_dir;

    ornt = vertices(iv).patch_reorientation(iptc,[1 2]);
    if (all (ornt == [0 0]))
      ind_sub{1} = 1:6;
      ind_sub{2} = 1:6;
    elseif (all (ornt == [1 0]))
      ind_sub{1} = ndof_dir_bsp_c(1)-5:ndof_dir_bsp_c(1);
      ind_sub{2} = 1:6;
    elseif (all (ornt == [0 1]))
      ind_sub{1} = 1:6;
      ind_sub{2} = ndof_dir_bsp_c(2)-5:ndof_dir_bsp_c(2);
    elseif (all (ornt == [1 1]))
      ind_sub{1} = ndof_dir_bsp_c(1)-5:ndof_dir_bsp_c(1);
      ind_sub{2} = ndof_dir_bsp_c(2)-5:ndof_dir_bsp_c(2);
    end
    for idim = 1:ndim
      [aux{idim},~] = find (Proj_patch{idim}(:,ind_sub{idim}));
      aux{idim} = unique (aux{idim});
    end
    [z{1:ndim}] = ndgrid (aux{:});
    auxI = sub2ind ([ndof_dir_bsp_f, 1], z{:});
    local_interior_dofs = sp_get_local_interior_functions (hspace.space_of_level(lev+1), patches(iptc));
    [~,local_indices,~] = intersect (local_interior_dofs, auxI);
%     [~,local_indices,~] = intersect (hspace.space_of_level(lev+1).interior_dofs_per_patch{patches(iptc)}, auxI);
    [~,interior_dofs_patch] = sp_get_functions_on_patch (hspace.space_of_level(lev+1), patches(iptc));
    possible_children = union (possible_children, interior_dofs_patch(local_indices(:)));
%     possible_children = union (possible_children, hspace.space_of_level(lev+1).dofs_on_patch{patches(iptc)}(local_indices(:)));
  end
  
% %  possible_children = cell2mat (hspace.space_of_level(lev+1).dofs_on_patch(patches));
  possible_children = union (possible_children, cell2mat (hspace.space_of_level(lev+1).dofs_on_edge(edges)));
  possible_children = union (possible_children, hspace.space_of_level(lev+1).dofs_on_vertex{iv});
  ref_matrix = matrix_basis_change__ (hspace, lev+1, ind_on_vertex, possible_children);

  for ii = 1:numel(ind_on_vertex)
    [auxI,~]= find(ref_matrix(:,ii));
    children_of_this_function = possible_children(auxI(:));
    children = union (children, children_of_this_function);
    children_of_function{position(ii)} = children_of_this_function(:).';  
  end
end

% % Children of edge and vertex functions. It uses a lot of memory
% [indices_not_interior, position] = setdiff (ind, 1:hspace.space_of_level(lev).ndof_interior);
% 
% ref_matrix = matrix_basis_change__ (hspace, lev+1, indices_not_interior);
% for ii = 1:numel(indices_not_interior)
%   [auxI,~]= find(ref_matrix(:,ind(position(ii))));
%   children = union (children, auxI);
%   children_of_function{position(ii)} = auxI;  
% end

if (nargout >= 2)
  flag = all (ismember (ind, hspace.active{lev}));
  if (~flag)
    flag = 2 *all (ismember (ind, union (hspace.active{lev}, hspace.deactivated{lev})));
  end
end

end
