function Y = my_kron(X)
%
% function Y = my_kron(X)
% 
% Sucessive application of kron.
%
% Input: X cell array, where X{idir} coordinates in the idir-th direction
%
% Output: Y column 
%
%

dim = numel(X);

switch dim
    case 1, Y = X{1}(:);
    case 2, Y = kron(X{2}(:),X{1}(:));
    case 3, Y = kron(X{3}(:),kron(X{2}(:),X{1}(:)));
end