function [h, meshsize] = get_meshsize(hmsh)

% Dominio parametrico


meshsize = zeros(hmsh.nlevels,1);
for lev = 1:hmsh.nlevels
ms = 1./hmsh.mesh_of_level(lev).nel_dir;
meshsize(lev) = norm(ms);
end

h = meshsize(hmsh.globnum_active(:,1));
