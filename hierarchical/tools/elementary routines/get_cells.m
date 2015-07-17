function I = get_cells(ind, degree, nelem)
%
% function I = get_cells(ind, degree, nelem)
%
% Computation of the indices of the cells within the support of a
% tensor-product B-spline function.
%
% Input:
%           ind: index of the B-spline basis function (row vector)
%           degree: polynomial degree in each coordinate direction
%           nelem: number of elements in each coordinate direction 
%
% Ouput:
%           I: matrix whose rows contains the indices of the cells in the
%           support of the given B-spline
%
% WARNING: We are assuming open knot vectors with multiplicity 1 at
% interior nodes (We generalize this function later...)
%
% Mejorar el siguiente procedimiento
%

if any(ind > (degree + nelem))
    I=[];
    disp('ERROR in function get_cells: Incorrect index for degree of freedom')
    return;
end    

dim = numel(ind);
y = cell(1,dim);
for i =1:dim
    y{i} = max(1,ind(i)-degree(i)):min(nelem(i),ind(i));
end
switch dim
    case 1, I = y{1}(:);
    case 2, [z{1:dim}] = ndgrid(y{1},y{2});
        I = [z{1}(:) z{2}(:)];
    case 3, [z{1:dim}] = ndgrid(y{1},y{2},y{3});
        I = [z{1}(:) z{2}(:) z{3}(:)];
end
