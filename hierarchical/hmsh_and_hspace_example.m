% prueba
clear all
close all

n = 6; % number of initial elements in each direction
degree = 2; % polynomial degree
dim = 2; % number of parametric directions

cells{1} = [3 3; 3 4; 4 3; 4 4];
[hmsh,hspace] = build_hspace_from_cells(dim,degree, n, cells, 0, true, true);

