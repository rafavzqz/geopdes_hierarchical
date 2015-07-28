% ex_hmsh_and_hspace
%
% Archivo de prueba, usa build_hspace_from_cells para construir una malla y
% un espacio jerarquicos a partir de las celdas desactivadas de cada nivel
% 

clear all
close all

n = 6; % number of initial elements in each direction
degree = 2; % polynomial degree
dim = 2; % number of parametric directions

%cells{1} = [3 3; 3 4; 4 3; 4 4];
cells{1} = [15 16 21 22]';

% [hmsh, hspace] = build_hspace_from_cells(dim, p, initial_num_el, cells,space_type, boundary, graficar_malla)
[hmsh,hspace] = build_hspace_from_cells(dim,degree, n, cells, 1, true, true);

