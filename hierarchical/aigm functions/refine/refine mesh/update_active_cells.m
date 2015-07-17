function [E, EE, new_cells] = update_active_cells(E, EE, cells_to_ref, indE)
%
% function [E, EE, new_cells] = update_active_cells(E, EE, cells_to_ref, indE)
%
% This function updates the active cells (removes cells of level lev and add their children)
%
% Input:    E: (matrix) global indices of active cells of a fixed level (one row per cell)
%           EE: analogous to E for level lev + 1
%           cells_to_ref: global indices of active cells to be refined
%           indE: indices array such that cells_to_ref = E(indE,:)
%           lev: level
% 
% Output:   [E, EE] updated
%            new_cells: Global indices of the new cells of level lev+1
%
%
% This function uses:   split_cell

dim = size(cells_to_ref,2);
% Refine cells
%%%%%%%%%%%%%%%%%%%%%%
% Mejorar el siguiente bloque
%%%%%%%%%%%%%%%%%%%%%%
n_refined_cells = size(cells_to_ref,1);
n_new_cells = (2^dim)*n_refined_cells;
naux = size(EE,1); % number of active elements at level lev + 1
EE = [EE; zeros(n_new_cells, dim)];

new_cells = zeros(0,dim);
naux = naux +1;
for i = 1: n_refined_cells
    new_ind = split_cell(cells_to_ref(i,:));
    new_cells = union(new_cells,new_ind,'rows');
    for child = 1:(2^dim)
        EE(naux,:) = new_ind(child,:);
        naux = naux +1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%

% Update E by removing the cells to be refined
E(indE,:) = [];
