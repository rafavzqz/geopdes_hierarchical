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

% Mejorar este chequeo
if size(hmsh.active{1},2)~=size(M{1},2)
    disp('Error: Bad call to update_active_cells');
    return,
end

% We update hsmh.active and hmsh.deactivated
nlevels = hmsh.nlevels;

if ~isempty(M{nlevels}) % if a new level is going to be activated
    hmsh.nlevels = hmsh.nlevels + 1;
    hmsh.active{nlevels+1,1} = zeros(0,1);
    hmsh.deactivated{nlevels+1,1} = zeros(0,1);
    new_cells = cell(nlevels+1,1);
else
    new_cells = cell(nlevels,1);
end

for i = 1:numel(new_cells)
    new_cells{i} = zeros(0,1);
end

for lev = 1:nlevels
    if ~isempty(M{lev})
        if isempty(indices)
            [unos, indE] = ismember(M{lev},hmsh.active{lev},'rows'); % M{lev} = hmsh.active{lev}(indE)
            if any(unos~=1)
                disp('Warning: update_active_cells: Some nonactive cells were selected');
            end
        else
            indE = indices{lev};
        end
        % Update hmsh.active{lev} by removing the cells to be refined
        hmsh.active{lev}(indE,:) = [];
        
        % Compute the children of cells in M{lev}
        % El siguiente loop se puede eliminar con una nueva version de
        % split_cell, que reciba varias celdas a la vez
        for i = 1: numel(M{lev})
            new_ind = split_cell(hmsh, lev, M{lev}(i));
            new_cells{lev+1} = union(new_cells{lev+1}, new_ind,'rows');
        end
        
        % Update hmsh.deactivated{lev} by adding the cells that were deactivated
        hmsh.deactivated{lev} = union (hmsh.deactivated{lev}, M{lev});
    end
end % for lev = 1:nlevels

for lev = 1:nlevels
    if ~isempty(M{lev})
        % Update hmsh.active{lev+1} by adding the children of the cells in M{lev},
        % i.e., the new cells of level lev+1
        hmsh.active{lev+1} = union(hmsh.active{lev+1}, new_cells{lev+1},'rows');
    end
end % for lev = 1:nlevels


% We update hmsh.nel_per_level
nel_per_level = zeros(1, hmsh.nlevels);
%aux_ind_el = [];
for lev = 1:hmsh.nlevels
    nel_per_level(lev) = size(hmsh.active{lev},1);
 %   aux_ind_el = vertcat(aux_ind_el, lev*ones(nel_per_level(lev),1));
end % for lev = 1:nlevels

hmsh.nel_per_level = nel_per_level;

% We update hmsh.nel
hmsh.nel = sum(hmsh.nel_per_level);

% We update hmsh.globnum_active
hmsh.globnum_active = zeros(hmsh.nel,hmsh.ndim+1);
Ne = cumsum([0; hmsh.nel_per_level(:)]);
for lev = 1:hmsh.nlevels
    ind_e = (Ne(lev)+1):Ne(lev+1);
    if ~isempty(ind_e)
        hmsh.globnum_active(ind_e, 1) = lev;
        globnum = cell(1,hmsh.ndim);
        [globnum{:}] = ind2sub(hmsh.mesh_of_level(lev).nel_dir,hmsh.active{lev}(:));
        hmsh.globnum_active(ind_e, 2:end) = cell2mat(globnum);
    end
end
end


function I = split_cell(hmsh, lev, ind)
%
% function I = split_cell(hmsh, lev, ind)
%
% Split a cell by dyadic refinement.
%
% Input:
%           ind: index of the cell
%
% Ouput:
%           I: column array that contains the indices of the children of the
%           given cell
%

% Mejorar el siguiente procedimiento
%


cells_sub = cell(1,hmsh.ndim);
[cells_sub{:}] = ind2sub(hmsh.mesh_of_level(lev).nel_dir, ind);
cells_sub = cell2mat(cells_sub);

% DYADIC REFINEMENT
aux = [2*cells_sub-1; 2*cells_sub];

switch hmsh.ndim
    case 1, I = aux;
    case 2, [z{1:hmsh.ndim}] = ndgrid(aux(:,1),aux(:,2));
        %I = [z{1}(:) z{2}(:)];
        I = sub2ind(hmsh.mesh_of_level(lev+1).nel_dir,z{1}(:), z{2}(:));
    case 3, [z{1:hmsh.ndim}] = ndgrid(aux(:,1),aux(:,2),aux(:,3));
        % I = [z{1}(:) z{2}(:) z{3}(:)];
        I = sub2ind(hmsh.mesh_of_level(lev+1).nel_dir,z{1}(:), z{2}(:),z{3}(:));
end

end
