function [hmsh, new_cells] = update_active_cells(hmsh, M, indices)
%
% function [hmsh, new_cells] = update_active_cells(hmsh, M, indices)
%
% This function updates the active cells (hmsh.active and hmsh.globnum_active) and deactivated cells (hmsh.deactivated) in each level when
% refining the cells in M. This function also updates hmsh.nlevels, hmsh.nel and hmsh.nel_per_level
%
% Input:    hmsh:
%           M{lev}: global indices of marked cells of level lev (one row per cell)
%           indices{lev}: array such that M{lev} =
%           hmsh.active{lev}(indices{lev}). If indices is [], then this function computes this information.
%
% Output:   hmsh:
%           new_cells{lev}: global indices of the new cells of level lev (one row per cell)
%

% This function uses split_cell, that is defined below
%

% We update hsmh.active and hmsh.deactivated
nlevels = hmsh.nlevels;

if (~isempty(M{nlevels})) % if a new level is going to be activated
  hmsh.nlevels = hmsh.nlevels + 1;
  hmsh.active{nlevels+1} = [];
  hmsh.deactivated{nlevels+1} = [];
  new_cells = cell (nlevels+1, 1);
else
  new_cells = cell (nlevels, 1);
end

for lev = 1:nlevels
  if (~isempty(M{lev}))
    if (isempty(indices))
      [unos, indE] = ismember (M{lev}, hmsh.active{lev});
      if (~all (unos))
        disp('Warning: update_active_cells: Some nonactive cells were selected');
      end
    else
      indE = indices{lev};
    end
    % Update hmsh.active{lev} by removing the cells to be refined
    hmsh.active{lev}(indE) = [];
      
    % Compute the children of cells in M{lev}
    % El siguiente loop se puede eliminar con una nueva version de
    % split_cell, que reciba varias celdas a la vez
    new_cells{lev+1} = split_cells (hmsh, lev, M{lev});
    % Update hmsh.deactivated{lev} by adding the cells that were deactivated
    hmsh.deactivated{lev} = union (hmsh.deactivated{lev}, M{lev});
  end
end

for lev = 1:nlevels
  if (~isempty(M{lev}))
    % Update hmsh.active{lev+1} by adding the children of the cells in M{lev},
    % i.e., the new cells of level lev+1
    hmsh.active{lev+1} = union (hmsh.active{lev+1}, new_cells{lev+1});
  end
end

% We update hmsh.nel_per_level and hmsh.nel
nel_per_level = zeros (1, hmsh.nlevels);
for lev = 1:hmsh.nlevels
  nel_per_level(lev) = numel (hmsh.active{lev});
end
hmsh.nel_per_level = nel_per_level;
hmsh.nel = sum (nel_per_level);

% We update hmsh.globnum_active
hmsh.globnum_active = zeros (hmsh.nel,hmsh.ndim+1);
Ne = cumsum ([0; hmsh.nel_per_level(:)]);
for lev = 1:hmsh.nlevels
  ind_e = (Ne(lev)+1):Ne(lev+1);
  if (~isempty(ind_e))
    hmsh.globnum_active(ind_e, 1) = lev;
    globnum = cell (1,hmsh.ndim);
    [globnum{:}] = ind2sub (hmsh.mesh_of_level(lev).nel_dir, hmsh.active{lev}(:));
    hmsh.globnum_active(ind_e, 2:end) = cell2mat (globnum);
  end
end

end




function I = split_cells (hmsh, lev, ind)
%
% function I = split_cells (hmsh, lev, ind)
%
% Split a set of cells of hmsh, of level lev, with the subdivision given by hmsh.nsub.
%
% Input:
%     hmsh: the hierarchical mesh
%     lev:  level of the cells to subdivide
%     ind:  indices of the cells in the Cartesian grid
%
% Output:
%           I: indices of the children, with the numbering of the Cartesian grid
%

z = cell (hmsh.ndim, 1);
cells_sub = cell (hmsh.ndim, 1);
[cells_sub{:}] = ind2sub ([hmsh.mesh_of_level(lev).nel_dir, 1], ind); % The extra 1 makes it work in any dimension

I = [];
for ii = 1:numel(cells_sub{1})
  aux = cell (hmsh.ndim, 1);
  for idim = 1:hmsh.ndim
    aux{idim} = hmsh.nsub(idim)*(cells_sub{idim}(ii)-1)+1:hmsh.nsub(idim)*(cells_sub{idim}(ii));
  end
  [z{1:hmsh.ndim}] = ndgrid (aux{:});
  I = union (I, sub2ind ([hmsh.mesh_of_level(lev+1).nel_dir, 1], z{:}));
end

end
