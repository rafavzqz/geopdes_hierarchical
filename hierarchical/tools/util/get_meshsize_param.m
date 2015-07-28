function [h, meshsize] = get_meshsize_param(hmsh)
%
% function [h, meshsize] = get_meshsize_param(hmsh)
%
% Computation of the meshsize in the parametric domain, assuming underlying uniform
% tensor product meshes
%
% INPUT:    hmsh
%
% OUTPUT:   h(el) = diameter of el, for el = 1:hmsh.nel
%           meshsize(lev) = diameter of the cells of level lev, for lev = 1:hmsh.nlevels
%
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% ATENCION: En esta funcion se asume que cada grilla producto tensor es una
% malla uniforme. 
%


meshsize = zeros(hmsh.nlevels,1);
for lev = 1:hmsh.nlevels
ms = 1./hmsh.mesh_of_level(lev).nel_dir;
meshsize(lev) = norm(ms);
end

h = meshsize(hmsh.globnum_active(:,1));
