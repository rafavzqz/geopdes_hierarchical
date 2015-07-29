function u = assemble_and_solve(hmsh, hspace, problem_data)

%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_diff:       diffusion coefficient (epsilon in the equation)
%    - f:            source term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition

tic
disp('Assembling and solving the system:')
stiffness = op_gradu_gradv_hier (hspace, hspace, hmsh, problem_data.c_diff);
rhs = op_f_v_hier (hspace, hmsh, problem_data.f);
mass = op_u_v_hier (hspace, hspace, hmsh, problem_data.c_diff);

% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX TO BE TESTED
% Apply Neumann boundary conditions 
for iside = problem_data.nmnn_sides
  if (hmsh.ndim > 1)
% Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
    gside = @(varargin) problem_data.g(varargin{:},iside);
    dofs = hspace.boundary(iside).dofs;
    rhs(dofs) = rhs(dofs) + op_f_v_hier (hspace.boundary(iside), hmsh.boundary(iside), gside);
  else
    if (iside == 1)
      x = hmsh.mesh_of_level(1).breaks{1}(1);
    else
      x = hmsh.mesh_of_level(1).breaks{1}(end);
    end
    sp_side = hspace.boundary(iside);
    rhs(sp_side.dofs) = rhs(sp_side.dofs) + problem_data.g(x,iside);
  end
end

% Apply Dirichlet boundary conditions

u = zeros (hspace.ndof, 1);

[u_dirichlet, dirichlet_dofs] = hsp_drchlt_l2_proj (hspace, hmsh, problem_data.h, problem_data.drchlt_sides);
u(dirichlet_dofs) = u_dirichlet;

% dirichlet_dofs = unique(cat(1,hspace.boundary(problem_data.drchlt_sides).dofs));

int_dofs = setdiff (1:hspace.ndof, dirichlet_dofs);

rhs(int_dofs) = rhs(int_dofs) - stiffness(int_dofs, dirichlet_dofs)*u(dirichlet_dofs);

% Solve the linear system
u(int_dofs) = stiffness(int_dofs, int_dofs) \ rhs(int_dofs);

n_int_dofs = numel(int_dofs);

tempo = toc;

fprintf('Total DOFs: %d - Used DOFs: %d (%f seconds)\n', hspace.ndof, n_int_dofs, tempo);