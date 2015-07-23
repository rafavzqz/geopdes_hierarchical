function [E, R, new_cells] = update_active_cells(E, R, M, indices)
%
% function [E, R, new_cells] = update_active_cells(E, R, M, indices)
%
% This function updates the active cells (E) and deactived cells (R) in each level when
% refining the cells in M
%
% Input:    E{lev}: (matrix) global indices of active cells of level lev (one row per cell)
%           R{lev}: (matrix) global indices of deactived cells of level lev (one row per cell)
%           M{lev}: (matrix) global indices of marked cells of level lev (one row per cell)
%           indices{lev}: indices array such that M{lev} =
%           E{lev}(indices{lev},:). If indices is [], then this function computes this information.
%
% Output:   E{lev}: (matrix) global indices of active cells of level lev (one row per cell)
%           R{lev}: (matrix) global indices of deactived cells of level lev (one row per cell)
%           new_cells{lev}: (matrix) global indices of the new cells of level lev (one row per cell)
%
%
% This function uses:   split_cell
%

if size(E{1},2)~=size(M{1},2)
    disp('Error: Bad call to update_active_cells');
    return,
end

dim = size(E{1},2);

nlevels = numel(E);

if ~isempty(M{nlevels}) % if a new level is going to be activated
    E{nlevels+1,1} = zeros(0,dim);
    R{nlevels+1,1} = zeros(0,dim);
    new_cells = cell(nlevels+1,1);
else
    new_cells = cell(nlevels,1);
end

for i = 1:numel(new_cells)
    new_cells{i} = zeros(0,dim);
end

for lev = 1:nlevels
    if ~isempty(M{lev})
        if isempty(indices)
            [unos, indE] = ismember(M{lev},E{lev},'rows'); % M{lev} = E{lev}(indE,:)
            if any(unos~=1)
                disp('Warning: update_active_cells: Some nonactive cells were selected');
            end
        else
            indE = indices{lev};
        end
        % Update E{lev} by removing the cells to be refined
        E{lev}(indE,:) = [];
        
        % Compute the children of cells in M{lev}
        for i = 1: size(M{lev},1)
            new_ind = split_cell(M{lev}(i,:));
            new_cells{lev+1} = union(new_cells{lev+1}, new_ind,'rows');
        end
        
        % Update M{lev} by adding the cells that were deactivated
        R{lev} = union(R{lev}, M{lev},'rows');
    end
end % for lev = 1:nlevels

for lev = 1:nlevels
    if ~isempty(M{lev})
        % Update E{lev+1} by adding the children of the cells in M{lev},
        % i.e., the new cells of level lev+1
        E{lev+1} = union(E{lev+1}, new_cells{lev+1},'rows');
    end
end % for lev = 1:nlevels