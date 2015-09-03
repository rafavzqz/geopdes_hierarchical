% This function uses:    compute_functions_to_deactivate
%                        update_active_functions
%                        compute_matrices_for_changing_basis
%

% HSPACE_REFINE: refine the hierarchical space, updating the fields of the object.
%
%   hspace = hspace_refine (hspace, hmsh, marked, flag, new_cells)
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

% XXXXXXXXXX Move the function into the class. Change the name

function hspace = refine_hierarchical_space (hspace, hmsh, M, flag, new_cells)

boundary = ~isempty (hspace.boundary);

% Computation of indices of functions of level lev that will become
% nonactive when removing the functions or elements in M{lev}
M = compute_functions_to_deactivate (hmsh, hspace, M, flag);

% Computation of a tensor product space if a new level is activated,
%  and the 1D projectors between the previous level and the new one.
if (numel(hspace.space_of_level) < hmsh.nlevels)
  msh_level = hmsh.mesh_of_level(hmsh.nlevels);
  degree = hspace.space_of_level(hmsh.nlevels-1).degree;
  knots = kntrefine (hspace.space_of_level(hmsh.nlevels-1).knots, hmsh.nsub-1, degree, degree-1);
  hspace.space_of_level(hmsh.nlevels) = sp_bspline (knots, degree, msh_level);
  coarse_space = hspace.space_of_level(hmsh.nlevels-1).constructor (msh_level);

  for idim = 1:hmsh.ndim
    knt_coarse = coarse_space.knots{idim};
    knt_fine = hspace.space_of_level(hmsh.nlevels).knots{idim};
    hspace.Proj{hmsh.nlevels-1,idim} = basiskntins (degree(idim), knt_coarse, knt_fine);
  end
end

% Update of active functions
if (boundary && hmsh.ndim > 1)
  new_elem = new_cells.interior;
else
  new_elem = new_cells;
end

hspace = update_active_functions (hspace, hmsh, new_elem, M);


% Update C, the matrices for changing basis
C = cell (hmsh.nlevels, 1);
C{1} = speye (hspace.space_of_level(1).ndof);
C{1} = C{1}(:,hspace.active{1});

for lev = 2:hmsh.nlevels
  I = speye (hspace.space_of_level(lev).ndof);
  aux = 1;
  for idim = 1:hmsh.ndim
    aux = kron (hspace.Proj{lev-1,idim}, aux);
  end
  C{lev} = [aux*C{lev-1}, I(:,hspace.active{lev})];
end
hspace.C = C;
clear C

% Update of sp_lev
% XXXX This can be improved to save computational time
hspace.sp_lev = cell (hmsh.nlevels,1);
for ilev = 1:hmsh.nlevels
  hspace.sp_lev{ilev} = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'gradient', true, 'hessian', true);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill the information for the boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (boundary && hmsh.ndim > 1)
    for iside = 1:2*hmsh.ndim
        %%    ind =[2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind = [2 2 1 1] in 2D;
        ind = setdiff (1:hmsh.ndim, ceil(iside/2));
        i = setdiff (1:hmsh.ndim, ind);
        if mod(iside,2) == 1
            boundary_ind = ones(hspace.nlevels,1);
        else
            boundary_ind = zeros(hspace.nlevels,1);
            for lev = 1:hspace.nlevels
                boundary_ind(lev) = hspace.space_of_level(lev).ndof_dir(i);
            end
        end
        M_boundary = cell(size(M));
        for lev = 1:numel(M)
            % if ~isempty(M{lev})
            M_sub = cell(1,hmsh.ndim);
            [M_sub{:}] = ind2sub(hspace.space_of_level(lev).ndof_dir,  M{lev}(:));
            M_sub = cell2mat(M_sub);
            M_boundary{lev} = M_sub(M_sub(:,i) == boundary_ind(lev),ind);
            % Mejorar lo siguiente
            switch hmsh.boundary(iside).ndim %XXXXXX CHECK IF THIS IS CORRECT
                case 2,
                    M_boundary{lev} = sub2ind(hspace.boundary(iside).space_of_level(lev).ndof_dir,M_boundary{lev}(:,1),M_boundary{lev}(:,2));
                case 3,
                    M_boundary{lev} = sub2ind(hspace.boundary(iside).space_of_level(lev).ndof_dir,M_boundary{lev}(:,1),M_boundary{lev}(:,2),M_boundary{lev}(:,3));
            end
            %end
        end
        hspace.boundary(iside) = refine_hierarchical_space(hspace.boundary(iside), hmsh.boundary(iside), ...
            M_boundary, 'functions', new_cells.boundary{iside});
        % Now, we fill hspace.boundary(iside).dofs
        globnum_active_boundary = [hspace.boundary(iside).globnum_active(:,1:i) boundary_ind(hspace.boundary(iside).globnum_active(:,1)) ...
            hspace.boundary(iside).globnum_active(:,(i+1):end)];
        [unos, hspace.boundary(iside).dofs] = ismember(globnum_active_boundary,hspace.globnum_active,'rows');
        if any(unos~=1)
            disp('Warning: Error when computing hspace.boundary().dofs')
            pause
        end
    end
else
    hspace.boundary = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%