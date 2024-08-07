% ADAPTIVITY_SOLVE_LAPLACE_MP_C1: assemble and solve the linear system for
%  Laplacian problem, using hierarchical spaces of C^1 multipatch splines.
%
% The function solves the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^n)
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% USAGE:
%
% u = adaptivity_solve_laplace_mp_C1 (hmsh, hspace, problem_data, method_data)
%
% INPUT:
%
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh_mp)
%   hspace: object representing the space of hierarchical splines (see hierarchical_space_mp_C1)
%   problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - weak_drchlt_sides: sides with Dirichlet boundary condition
%    - c_diff:       diffusion coefficient (epsilon in the equation)
%    - f:            function handle of the source term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
%   method_data: a structure with data for discretization. For this function, it must contain the fields:
%    - Cpen: penalty parameter for Nitsche's method.
%
% OUTPUT:
%
%   u: computed degrees of freedom
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2022-2023 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function u = adaptivity_solve_laplace_mp_C1 (hmsh, hspace, problem_data, method_data)

data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end

% Compute and assemble the matrices 
stiff_mat = op_gradu_gradv_hier (hspace, hspace, hmsh, c_diff);
rhs = op_f_v_hier (hspace, hmsh, f);

% Apply Neumann boundary conditions
ndofs=0;
for ilev = 1:hmsh.nlevels
    ndofs = ndofs + hspace.ndof_per_level(ilev);
    for iref = nmnn_sides
        gref = @(varargin) g(varargin{:}, iref);
        for bnd_side = 1:hmsh.mesh_of_level(ilev).boundaries(iref).nsides
            iptc = hmsh.mesh_of_level(ilev).boundaries(iref).patches(bnd_side);
            iside = hmsh.mesh_of_level(ilev).boundaries(iref).faces(bnd_side);

            msh_side = hmsh.mesh_of_level(ilev).msh_patch{iptc}.boundary(iside);
            sp_side = hspace.space_of_level(ilev).sp_patch{iptc}.boundary(iside);
            rhs_nmnn = op_f_v_tp (sp_side, msh_side, gref);
            
            dofs = 1:ndofs;
            
            rhs_aux=hspace.space_of_level(ilev).Cpatch{iptc}(sp_side.dofs,:).' * rhs_nmnn; %all the dofs of the current level are included here (rows) 
            rhs(dofs) = rhs(dofs) + hspace.Csub{ilev}.' * rhs_aux;
        end
    end
end

% Apply Dirichlet boundary conditions in weak form, by Nitsche's method
if (exist ('weak_drchlt_sides', 'var'))
  if (isfield (method_data, 'Cpen'))
    [N_mat, N_rhs] = sp_weak_drchlt_bc_laplace (hspace, hmsh, weak_drchlt_sides, h, c_diff, method_data.Cpen); 
  else
    [N_mat, N_rhs] = sp_weak_drchlt_bc_laplace (hspace, hmsh, weak_drchlt_sides, h, c_diff);
  end
        
  stiff_mat = stiff_mat - N_mat;
  rhs = rhs + N_rhs;
end

% Solve the linear system
u = stiff_mat \ rhs;

end