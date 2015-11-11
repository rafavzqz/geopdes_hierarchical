% This function uses:    compute_functions_to_deactivate
%                        update_active_functions

% HSPACE_REFINE: refine the hierarchical space, updating the fields of the object.
%
%   [hspace, Cref] = hspace_refine (hspace, hmsh, marked, flag, new_cells)
%
% INPUT:
%
%   hspace:    object representing the coarse hierarchical space (see hierarchical_space)
%   hmsh:      object representing the refined hierarchical mesh (see hierarchical_mesh)
%   marked:    cell array with the indices, in the tensor product setting, of the marked elements/functions for each level
%   flag:      the refinement strategy, marking either 'elements' or 'functions'.
%   new_cells: cell array with the global indices of the new active elements for each level
%
% OUTPUT:
%
%   hspace:    object representing the refined hierarchical space (see hierarchical_space)
%   Cref:      a matrix to pass from the coarse space (input) to the refined space (output)
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
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

function [hspace,Cref] = hspace_refine (hspace, hmsh, M, flag, new_cells)

boundary = ~isempty (hspace.boundary);

% Computation of indices of functions of level lev that will become
% nonactive when removing the functions or elements in M{lev}
M = compute_functions_to_deactivate (hmsh, hspace, M, flag);

% Computation of a tensor product space if a new level is activated,
%  and the 1D projectors between the previous level and the new one.
if (numel(hspace.space_of_level) < hmsh.nlevels)
  msh_level = hmsh.mesh_of_level(hmsh.nlevels);
  
  if (isa (hspace.space_of_level(1), 'sp_scalar'))
    degree = hspace.space_of_level(hmsh.nlevels-1).degree;
    knots = kntrefine (hspace.space_of_level(hmsh.nlevels-1).knots, hmsh.nsub-1, degree, degree-1);
    hspace.space_of_level(hmsh.nlevels) = sp_bspline (knots, degree, msh_level);

    for idim = 1:hmsh.ndim
      knt_coarse = hspace.space_of_level(hmsh.nlevels-1).knots{idim};
      knt_fine = hspace.space_of_level(hmsh.nlevels).knots{idim};
      hspace.Proj{idim,hmsh.nlevels-1} = basiskntins (degree(idim), knt_coarse, knt_fine);
    end
    
  elseif (isa (hspace.space_of_level(1), 'sp_vector'))
    scalar_spaces = cell (hspace.ncomp, 1);
    for icomp = 1:hspace.ncomp
      sp_coarse = hspace.space_of_level(hmsh.nlevels-1).scalar_spaces{icomp};
      degree = sp_coarse.degree;
      knots = kntrefine (sp_coarse.knots, hmsh.nsub-1, degree, degree-1);
      scalar_spaces{icomp} = sp_bspline (knots, degree, msh_level);
      
      for idim = 1:hmsh.ndim
        knt_coarse = sp_coarse.knots{idim};
        knt_fine = knots{idim};
        hspace.Proj{icomp,idim,hmsh.nlevels-1} = basiskntins (degree(idim), knt_coarse, knt_fine);
      end
    end
    hspace.space_of_level(hmsh.nlevels) = sp_vector (scalar_spaces, msh_level);
  end

  hspace.nlevels = hmsh.nlevels;
  hspace.active{hmsh.nlevels} = [];
  hspace.deactivated{hmsh.nlevels} = [];
  hspace.ndof_per_level(hmsh.nlevels) = 0;
  
  M{hmsh.nlevels} = [];
end

% Update of active functions
[hspace,Cref] = update_active_functions (hspace, hmsh, new_cells, M);

% Update C, the matrices for changing basis
C = cell (hmsh.nlevels, 1);
C{1} = speye (hspace.space_of_level(1).ndof);
C{1} = C{1}(:,hspace.active{1});

for lev = 2:hmsh.nlevels
  I = speye (hspace.space_of_level(lev).ndof);
  aux = matrix_basis_change (hspace, hmsh, lev);
  C{lev} = [aux*C{lev-1}, I(:,hspace.active{lev})];
end
hspace.C = C;
clear C

% Fill the information for the boundaries
if (boundary)% && hmsh.ndim > 1)

  Nf = cumsum ([0, hspace.ndof_per_level]);
  for iside = 1:2*hmsh.ndim
    %%    ind =[2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind = [2 2 1 1] in 2D;
    %%    ind2 = [1 1 2 2 3 3] in 3D, ind2 = [1 1 2 2] in 2D
    ind2 = ceil (iside/2);
    if (mod(iside,2) == 1)
      boundary_ind = ones (hspace.nlevels,1);
    else
      boundary_ind = zeros (hspace.nlevels,1);
      for lev = 1:hspace.nlevels
        boundary_ind(lev) = hspace.space_of_level(lev).ndof_dir(ind2);
      end
    end

    if (hmsh.ndim > 1)
      M_boundary = cell (size (M));
      for lev = 1:numel (M)
        [~,~,M_boundary{lev}] = intersect (M{lev}, hspace.space_of_level(lev).boundary(iside).dofs);
      end
      
      new_cells_boundary = cell (size (new_cells));
      for lev = 1:numel (new_cells)
        new_cells_boundary{lev} = get_boundary_indices (iside, hmsh.mesh_of_level(lev).nel_dir, new_cells{lev});
      end
      hspace.boundary(iside) = hspace_refine (hspace.boundary(iside), hmsh.boundary(iside), ...
          M_boundary, 'functions', new_cells_boundary);

      dofs = [];
      for lev = 1:hspace.boundary(iside).nlevels
        [~,iact] = intersect (hspace.active{lev}, hspace.space_of_level(lev).boundary(iside).dofs);
        dofs = union (dofs, Nf(lev) + iact);
      end
      hspace.boundary(iside).dofs = dofs;

    elseif (hmsh.ndim == 1)
      dofs = [];
      for lev = 1:hspace.nlevels
        [~,iact] = intersect (hspace.active{lev}, hspace.space_of_level(lev).boundary(iside).dofs);
        dofs = union (dofs, Nf(lev) + iact);
      end
      hspace.boundary(iside).dofs = dofs;
    end
  end
  
else
    hspace.boundary = [];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function hspace = update_active_functions(hspace, hmsh, new_cells, marked_fun)
%
% This function updates the active dofs (hspace.active), their coefficients (hspace.coeff_pou) and deactivated dofs (hspace.deactivated) in each level when
% refining the functions in marked_fun. This function also updates hspace.nlevels, hspace.ndof and hspace.ndof_per_level
%
% Input:    hspace:    the coarse mesh, an object of the class hierarchical_space
%           hmsh:      an object of the class hierarchical_mesh, already refined
%           new_cells: cells added to the refined mesh, see hmsh_refine
%           marked_fun{lev}: indices of active functions of level lev to be deactivated
%
% Output:   hspace:    the refined mesh, an object of the class hierarchical_space
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
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

function [hspace,Cref] = update_active_functions (hspace, hmsh, new_cells, marked_fun)

active = hspace.active;
deactivated = hspace.deactivated;


for lev = 1:hspace.nlevels-1

% Remove the marked functions from the active functions of level lev
  active{lev} = setdiff (active{lev}, marked_fun{lev});
  deactivated{lev} = union (marked_fun{lev}, deactivated{lev});

  if (strcmpi (hspace.type, 'simplified') && ~isempty (marked_fun{lev}))

    Cmat = matrix_basis_change (hspace, hmsh, lev+1);  
    [ii,~] = find (Cmat(:,marked_fun{lev}));

    active_and_deact = union (active{lev+1}, deactivated{lev+1});
    new_active = setdiff (unique (ii), active_and_deact);
    active{lev+1} = union (active{lev+1}, new_active);

% Mark functions whose support has been already refined completely
    [~, cells_per_fun] = sp_get_cells (hspace.space_of_level(lev+1), hmsh.mesh_of_level(lev+1), new_active);
    flag = cellfun (@(x) isempty (intersect (x, hmsh.active{lev+1})), cells_per_fun);
    marked_fun{lev+1} = union (marked_fun{lev+1}, new_active(flag==1,:));

  elseif (strcmpi (hspace.type, 'standard') && ~isempty (new_cells{lev+1}))
    new_possible_active_fun = sp_get_basis_functions (hspace.space_of_level(lev+1), hmsh.mesh_of_level(lev+1), new_cells{lev+1});
    new_possible_active_fun = setdiff (new_possible_active_fun, active{lev+1});

    [~, elem] = sp_get_cells (hspace.space_of_level(lev+1), hmsh.mesh_of_level(lev+1), new_possible_active_fun);

    new_functions = cellfun (@(x) all (ismember (x, union (hmsh.active{lev+1}, hmsh.deactivated{lev+1}))), elem);
    active{lev+1} = union (active{lev+1}, new_possible_active_fun(new_functions));
  end
  
end % for lev


%XXXXXXXXXXX The computation of Cref for truncated functions is wrong!!
% Computation of the matrix to pass from the original to the refined space
if (nargout == 2 || ~hspace.truncated)
  ndlev = hspace.ndof_per_level(1);
  active_and_deact = union (active{1}, deactivated{1});
  [~,indices] = intersect (active_and_deact, hspace.active{1});
  Id = sparse (numel(active_and_deact), ndlev);
  Id(indices,:) = speye (ndlev, ndlev);
  Cref = Id;

  for lev = 1:hspace.nlevels-1

    Cmat = matrix_basis_change (hspace, hmsh, lev+1);

    [~,deact_indices] = intersect (active_and_deact, deactivated{lev});
  
    aux = Cref;

    ndof_per_level = cellfun (@numel, active);
    ndof_prev_levs = sum (ndof_per_level(1:lev-1));
    [~,~,indices] = intersect (active{lev}, active_and_deact);
    Cref(ndof_prev_levs+(1:numel(active{lev})),:) = Cref(ndof_prev_levs+indices,:); %ignoring deactivated{lev}

    active_and_deact = union (active{lev+1}, deactivated{lev+1});
  
    ndof_until_lev = sum (ndof_per_level(1:lev));
    Cref(ndof_until_lev+(1:numel(active_and_deact)),:) = ...
      Cmat(active_and_deact,deactivated{lev}) * aux(ndof_prev_levs+deact_indices,:); %deactivated{lev} in terms of active_deact{lev+1}

    ndlev = hspace.ndof_per_level(lev+1);
    [~,indices] = intersect (active_and_deact, hspace.active{lev+1});
    Id = sparse (numel(active_and_deact), ndlev);
    Id(indices,:) = speye (ndlev, ndlev);
    Cref = [Cref, [sparse(ndof_until_lev,ndlev); Id]];  
  
  end % for lev
end


hspace.active = active(1:hspace.nlevels);
hspace.deactivated = deactivated(1:hspace.nlevels);
hspace.ndof_per_level = cellfun (@numel, hspace.active);
hspace.ndof = sum (hspace.ndof_per_level);

if (hspace.truncated)
  hspace.coeff_pou = ones (hspace.ndof, 1);
else
  hspace.coeff_pou = Cref * hspace.coeff_pou;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function C = matrix_basis_change (hspace, hmsh, lev)
%
% Compute the new matrices to represent functions of the previous level 
% as linear combinations of splines (active and inactive) of the current level 
%
% Input:  hspace: an object of the class hierarchical_space
%         hmsh:   an object of the class hierarchical_mesh
%         lev:    the level for which we compute the matrix
%
% Output:   C:    matrix to change basis from level lev-1 to level lev

function C = matrix_basis_change (hspace, hmsh, lev)

if (isa (hspace.space_of_level(1), 'sp_scalar'))
  C = 1;
  for idim = 1:hmsh.ndim
    C = kron (hspace.Proj{idim,lev-1}, C);
  end
elseif (isa (hspace.space_of_level(1), 'sp_vector'))
  for icomp = 1:hspace.ncomp
    Caux{icomp} = 1;
    for idim = 1:hmsh.ndim
      Caux{icomp} = kron (hspace.Proj{icomp,idim,lev-1}, Caux{icomp});
    end
    C = blkdiag (Caux{:});
  end
end

if (hspace.truncated)
  indices = union (hspace.active{lev}, hspace.deactivated{lev});
  C(indices,:) = 0;
end

end
