function plot_cells(indices, hmsh, color)
%
% function plot_cells(indices, hmsh, color)
%
% Funcion personal --- Solo 2d
%

if nargin == 2
    color = 'b';
end

plot_msh(hmsh);
hold on
ncells = size (indices,1);
for i = 1:ncells
    lev = indices(i,1);
    mx = hmsh.mesh_of_level(lev).nel_dir(1);
    indx = indices(i,2);
    x = [(indx-1)/mx indx/mx];
    my = hmsh.mesh_of_level(lev).nel_dir(2);
    indy = indices(i,3);
    y = [(indy-1)/my indy/my];
    abcisas = [x(1) x(2) x(2) x(1) x(1)];
    ordenadas = [y(1) y(1) y(2) y(2) y(1)];
    fill(abcisas,ordenadas,color);
end
hold off
