function I = split_cell(ind)
%
% function I = split_cell(ind)
%
% Split a cell by dyadic refinement.
%
% Input:
%           ind: index of the cell (row vector)
%
% Ouput:
%           I: matrix whose rows contains the indices of the children of the
%           given cell
%
%
% Mejorar el siguiente procedimiento
%

dim = numel(ind);
aux = [2*ind-1; 2*ind];

switch dim
    case 1, I = aux;
    case 2, [z{1:dim}] = ndgrid(aux(:,1),aux(:,2));
        I = [z{1}(:) z{2}(:)];
    case 3, [z{1:dim}] = ndgrid(aux(:,1),aux(:,2),aux(:,3));
        I = [z{1}(:) z{2}(:) z{3}(:)];
end
