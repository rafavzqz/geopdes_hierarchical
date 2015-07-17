function [u_dirichlet, dirichlet_dofs] = hsp_drchlt_l2_proj (hspace, hmsh, h, drchlt_sides)

u_dir = [];
dir_dofs = [];

for iside = drchlt_sides
% Restrict the function handle to the specified side, in any dimension, hside = @(x,y) h(x,y,iside)
    hside = @(varargin) h(varargin{:},iside);
[rhs,stiffness, mass] = assemble (hmsh.boundary(iside), hspace.boundary(iside), hside);
u_side = mass\rhs;
dofs = hspace.boundary(iside).dofs;

u_dir = [u_dir; u_side];
dir_dofs = [dir_dofs; dofs]; 

end

[dirichlet_dofs, ind] = unique(dir_dofs);
u_dirichlet = u_dir(ind);

% Deberiamos verificar si u coincide en la interseccion de bordes!

