function Y = my_ndgrid(X)
%
% function Y = my_ndgrid(X)
% 
% Input: X cell array, where X{idir} coordinates in the idir-th direction
%
% Output: Y matrix, where each row contains the dim coordinates of a point 
%
% En realidad esta funcion se puede reemplazar usando ind2sub (cuando se usa para indices enteros)...
% Hasta ahora, solo se uso en main.m para inicializar variables globales
%
% Tambien se usa para armar la tabla de puntos de cuadratura en 2d o 3d, cuando se conocen los valores en cada dirección
%

dim = numel(X);

switch dim
    case 1, Y = X{1}(:);
    case 2, [Y{1:dim}] = ndgrid(X{1},X{2});
        Y = [Y{1}(:) Y{2}(:)];
    case 3, [Y{1:dim}] = ndgrid(X{1},X{2},X{3});
        Y = [Y{1}(:) Y{2}(:) Y{3}(:)];
end
