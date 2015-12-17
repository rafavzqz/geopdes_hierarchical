% HSPACE_REFINE: refine the hierarchical space, updating the fields of the object.
%
%   [hspace, Cref] = hspace_refine (hspace, hmsh, marked, flag, new_cells)
%
% INPUT:
%
%   hspace:    object representing the coarse hierarchical space (see hierarchical_space_mp)
%   hmsh:      object representing the refined hierarchical mesh (see hierarchical_mesh_mp)
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
  degree = hspace.space_of_level(hmsh.nlevels-1).sp_patch{1}.degree;
  [new_space, Proj] = sp_refine (hspace.space_of_level(hmsh.nlevels-1), msh_level, hmsh.nsub, degree, degree-1);
  hspace.space_of_level(hmsh.nlevels) = new_space; clear new_space
  hspace.Proj(hmsh.nlevels-1,:) = Proj(:);

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
  if (hmsh.ndim > 1)
    new_cells_boundary = cell (size (new_cells));
    levels = find (~cellfun (@isempty, new_cells));
    for lev = levels(:).'
      Nelem = cumsum ([0 hmsh.mesh_of_level(lev).nel_per_patch]);
      Nelem_bnd = cumsum ([0 hmsh.boundary.mesh_of_level(lev).nel_per_patch]);
      for iptc = 1:hmsh.boundary.npatch
        patch_number = hmsh.boundary.mesh_of_level(1).patch_numbers(iptc);
        side_number  = hmsh.boundary.mesh_of_level(1).side_numbers(iptc);
        [~,indices,~] = intersect (Nelem(patch_number)+1:Nelem(patch_number+1), new_cells{lev});
        bnd_indices = get_boundary_indices (side_number, hmsh.mesh_of_level(lev).msh_patch{patch_number}.nel_dir, indices);
        new_cells_boundary{lev} = union (new_cells_boundary{lev}, bnd_indices+Nelem_bnd(iptc));
      end
    end
    M_boundary = cell (size (M));
    levels = find (~cellfun (@isempty, M));
    for lev = levels(:).'
      [~,~,M_boundary{lev}] = intersect (M{lev}, hspace.space_of_level(lev).boundary.dofs);
    end
    
    hspace.boundary = hspace_refine (hspace.boundary, hmsh.boundary, M_boundary, 'functions', new_cells_boundary);

    Nf = cumsum ([0, hspace.ndof_per_level]);
    dofs = [];
    for lev = 1:hspace.boundary.nlevels
      [~,iact,jact] = intersect (hspace.active{lev}, hspace.space_of_level(lev).boundary.dofs);
      [~, reorder] = sort (jact);
      dofs = vertcat (dofs, Nf(lev) + iact(reorder));
    end
    hspace.boundary.dofs = dofs;
  elseif (hmsh.ndim == 1)
    error ('The 1D multipatch has not been implemented')
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
% Input:    hspace:    the coarse mesh, an object of the class hierarchical_space_mp
%           hmsh:      an object of the class hierarchical_mesh_mp, already refined
%           new_cells: cells added to the refined mesh, see hmsh_refine
%           marked_fun{lev}: indices of active functions of level lev to be deactivated
%
% Output:   hspace:    the refined mesh, an object of the class hierarchical_space_mp
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


% Computation of the matrix to pass from the original to the refined space
if (nargout == 2 || ~hspace.truncated)
    
  %THB case
  if (hspace.truncated)
      
  ndlev = hspace.ndof_per_level(1);
  tp_ndlev=hspace.space_of_level(1).ndof;
  active_and_deact = union (active{1}, deactivated{1});
  [~,indices] = intersect (active_and_deact, hspace.active{1});
  Id = sparse (tp_ndlev, ndlev);
  Id(indices,:) = speye (ndlev, ndlev);
  Cref = Id;

  for lev = 1:hspace.nlevels-1
    
    Cmat = matrix_basis_change(hspace, hmsh, lev+1); 
    n_act_deact=numel(active_and_deact);
  
    aux = Cref;

    ndof_per_level = cellfun (@numel, active);
    ndof_prev_levs = sum (ndof_per_level(1:lev-1));
    Cref(ndof_prev_levs+(1:numel(active{lev})),:) = Cref(ndof_prev_levs+active{lev},:);

    active_and_deact = union (active{lev+1}, deactivated{lev+1});
    tp_ndlev=hspace.space_of_level(lev+1).ndof;
  
    ndof_until_lev = sum (ndof_per_level(1:lev));
    Cref(ndof_until_lev+(1:tp_ndlev),:) = Cmat*aux(ndof_prev_levs+1:end,:);

    ndlev = hspace.ndof_per_level(lev+1);
    Id = sparse (tp_ndlev, ndlev);
    Id(hspace.active{lev+1},:) = speye (ndlev, ndlev);
    Cref = [Cref, [sparse(ndof_until_lev,ndlev); Id]];  
   
  end
  ndof_per_level = cellfun (@numel, active);
  ndof_prev_levs = sum (ndof_per_level(1:hspace.nlevels-1));
  [~,~,indices] = intersect (active{hspace.nlevels}, active_and_deact);
  Cref(ndof_prev_levs+(1:numel(active{hspace.nlevels})),:) = Cref(ndof_prev_levs+indices,:);
  Cref(ndof_prev_levs+numel(active{hspace.nlevels})+1:end,:)=[];
  
  else  
    
  %HB-splines case  
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
    Cref(ndof_prev_levs+(1:numel(active{lev})),:) = Cref(ndof_prev_levs+indices,:); 

    active_and_deact = union (active{lev+1}, deactivated{lev+1});
  
    ndof_until_lev = sum (ndof_per_level(1:lev));
    Cref(ndof_until_lev+(1:numel(active_and_deact)),:) = ...
      Cmat(active_and_deact,deactivated{lev}) * aux(ndof_prev_levs+deact_indices,:); 

    ndlev = hspace.ndof_per_level(lev+1);
    [~,indices] = intersect (active_and_deact, hspace.active{lev+1});
    Id = sparse (numel(active_and_deact), ndlev);
    Id(indices,:) = speye (ndlev, ndlev);
    Cref = [Cref, [sparse(ndof_until_lev,ndlev); Id]];  
   
  end  
  end
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
% Input:  hspace: an object of the class hierarchical_space_mp
%         hmsh:   an object of the class hierarchical_mesh_mp
%         lev:    the level for which we compute the matrix
%
% Output:   C:    matrix to change basis from level lev-1 to level lev

function C = matrix_basis_change (hspace, hmsh, lev)

C = sparse (hspace.space_of_level(lev).ndof, hspace.space_of_level(lev-1).ndof);
for ipatch = 1:hmsh.npatch
  Cpatch = 1;
  Proj = hspace.Proj{lev-1, ipatch};
  for idim = 1:hmsh.ndim
    Cpatch = kron (Proj{idim}, Cpatch);
  end
  
  if (strcmpi (hspace.space_of_level(1).sp_patch{ipatch}.space_type, 'NURBS'))
    sp_patch = hspace.space_of_level(lev-1).sp_patch{ipatch};
    Wlev = spdiags (sp_patch.weights(:), 0, sp_patch.ndof, sp_patch.ndof);
    sp_patch = hspace.space_of_level(lev).sp_patch{ipatch};
    Wlev_fine = spdiags (sp_patch.weights(:), 0, sp_patch.ndof, sp_patch.ndof);
    Cpatch = Wlev_fine * Cpatch * Wlev;
  end
  C(hspace.space_of_level(lev).gnum{ipatch}, hspace.space_of_level(lev-1).gnum{ipatch}) = Cpatch;
end


if (hspace.truncated)
  indices = union (hspace.active{lev}, hspace.deactivated{lev});
  C(indices,:) = 0;
end

end
