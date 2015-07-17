function I = get_basis_functions(ind, degree, nelem)
%
% function I = get_basis_functions(ind, degree, nelem)
%
% Computation of indices of the tensor-product B-splines acting on a cell.
%
% Input:
%           ind: index of the cell (row vector)
%           degree: polynomial degree in each coordinate direction
%           nelem: number of elements in each coordinate direction 
%
% Ouput:
%           I: matrix whose rows contains the indices of the B-spline basis functions acting on the
%           given cell
%
% WARNING: We are assuming open knot vectors with multiplicity 1 at
% interior nodes (We generalize this function later...)
%
% Mejorar el siguiente procedimiento
%

if any(ind > nelem)
    I=[];
    disp('ERROR in function get_basis_functions: Incorrect index for a cell')
    return;
end    

dim = numel(ind);
y = cell(1,dim);
for i =1:dim
    y{i} = ind(i):(ind(i)+degree(i));
end
switch dim
    case 1, I = y{1}(:);
    case 2, [z{1:dim}] = ndgrid(y{1},y{2});
        I = [z{1}(:) z{2}(:)];
    case 3, [z{1:dim}] = ndgrid(y{1},y{2},y{3});
        I = [z{1}(:) z{2}(:) z{3}(:)];
end
