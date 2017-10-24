% adaptivity_solve_nonlinear_generalized_alphaMethod: assemble and solve the linear system for
% Poisson transient problem, using hierarchical spaces and generlized alpha-method
% time integartion scheme
%
% The function solves the non-linear problem
%
%  c(x,u) d(u)/dt  - div ( epsilon(x,u) grad (u)) = f(t)    in Omega = F((0,1)^n)
%                              epsilon(x,u) du/dn = g       on Gamma_N
%                                               u = h       on Gamma_D
%
% USAGE:
%
% u = adaptivity_solve_nonlinear_generalized_alphaMethod (hmsh, hspace, method_data)
%
% INPUT:
%
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the space of hierarchical splines (see hierarchical_space)
%   u_prev: solution vector of the last iteration
%   problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_diff:       temperature dependent diffusion coefficient (epsilon in the equation)
%    - c_cap:        temperature dependent heat capacity (c in the equation)
%    - f:            function handle of the source term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
%
% OUTPUT:
%
%   u: computed degrees of freedom
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
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

function [u, problem_data] = adaptivity_solve_nonlinear_generalized_alphaMethod (hmsh, hspace, time_step, problem_data, plot_data, u_0)

iter = 0;

u = u_0;
res = 1.0;
problem_data.non_linear_convergence_flag = 0;
% assume homogeneus time discretization
delta_t = problem_data.time_discretization(2) - problem_data.time_discretization(1);

while (norm(res) > problem_data.Newton_tol)
    iter = iter + 1;
    
    if (plot_data.print_info)
        fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Newton-Raphson iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',iter);
    end
    
    stiff_mat = op_gradu_gradv_nl_hier (hspace, hspace, hmsh, problem_data.c_diff, u);
    mass_mat = op_u_v_nl_hier(hspace, hspace, hmsh, problem_data.c_cap, u, u_0);
    %lumped mass matrix
    if problem_data.lumped
        warning('Lumping the Mass Matrix might decrease the approximation quality !!!');
        mass_mat = diag(sum(mass_mat));
    end
    f = op_f_v_time_hier (hspace, hmsh, problem_data, time_step);
    if time_step > 1
        f_prev = op_f_v_time_hier (hspace, hmsh, problem_data, time_step-1);
    else
        f_prev = f;
    end
    
    % Apply Neumann boundary conditions
    if (~isfield (struct (hmsh), 'npatch')) % Single patch case
        for iside = problem_data.nmnn_sides
            if (hmsh.ndim > 1)
                % Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
                gside = @(varargin) problem_data.g(varargin{:},iside);
                dofs = hspace.boundary(iside).dofs;
                f(dofs) = f(dofs) + op_f_v_hier (hspace.boundary(iside), hmsh.boundary(iside), gside);
            else
                if (iside == 1)
                    x = hmsh.mesh_of_level(1).breaks{1}(1);
                else
                    x = hmsh.mesh_of_level(1).breaks{1}(end);
                end
                sp_side = hspace.boundary(iside);
                f(sp_side.dofs) = f(sp_side.dofs) + problem_data.g(x,iside);
            end
        end
    else % Multipatch case
        boundaries = hmsh.mesh_of_level(1).boundaries;
        Nbnd = cumsum ([0, boundaries.nsides]);
        for iref = problem_data.nmnn_sides
            iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
            gref = @(varargin) problem_data.g(varargin{:},iref);
            rhs_nmnn = op_f_v_hier (hspace.boundary, hmsh.boundary, gref, iref_patch_list);
            f(hspace.boundary.dofs) = f(hspace.boundary.dofs) + rhs_nmnn;
        end
    end
    
    % Apply convection boundary conditions
    for iside = problem_data.convection_sides
        % Restrict the function handle to the specified side, in any dimension, cside = @(x,y,z) conv(x,y,z,iside)
        cside = @(varargin) problem_data.conv_fun(varargin{:},iside);
        cside_rhs = @(varargin) problem_data.conv_fun_rhs(varargin{:},iside);
        dofs = hspace.boundary(iside).dofs;
        mass_mat(dofs,dofs) = mass_mat(dofs,dofs) + ...
            op_u_v_hier(hspace.boundary(iside), hspace.boundary(iside), hmsh.boundary(iside), cside);
        f(dofs) = f(dofs) + op_conv_v_hier (hspace.boundary(iside), hmsh.boundary(iside), cside_rhs, u);
    end
    
    % Apply radiation boundary conditions
    for iside = problem_data.radiation_sides
        % Restrict the function handle to the specified side, in any dimension, rside = @(x,y,z) rad(x,y,z,iside)
        rtside = @(varargin) problem_data.rad_tilda_fun(varargin{:},iside);
        rside_rhs = @(varargin) problem_data.rad_fun(varargin{:},iside);
        
        dofs = hspace.boundary(iside).dofs;
        mass_mat(dofs,dofs) = mass_mat(dofs,dofs) + ...
            op_u_v_hier(hspace.boundary(iside), hspace.boundary(iside), hmsh.boundary(iside), rtside);
        f(dofs) = f(dofs) + op_rad_v_hier (hspace.boundary(iside), hmsh.boundary(iside), rside_rhs, u);
    end
    
    % Apply Dirichlet boundary conditions
    if ~isempty(problem_data.drchlt_sides)
        [u_dirichlet, dirichlet_dofs] = sp_drchlt_l2_proj (hspace, hmsh, problem_data.h, problem_data.drchlt_sides);
        u(dirichlet_dofs) = u_dirichlet;
    else
        dirichlet_dofs = [];
    end
    
    int_dofs = setdiff (1:hspace.ndof, dirichlet_dofs);
    alpha = problem_data.alpha;
    %residuum
    if (plot_data.print_info)
        fprintf('\n \t Assembly residuum vector');
    end
    if ~isempty(dirichlet_dofs)
        f(int_dofs) = (alpha*f(int_dofs) - stiff_mat(int_dofs, dirichlet_dofs)*u(dirichlet_dofs) + (1-alpha)*f_prev(int_dofs))* delta_t +...
            mass_mat(int_dofs, int_dofs)*u_0(int_dofs) - mass_mat(int_dofs, dirichlet_dofs)*u(dirichlet_dofs);
    else
        f(int_dofs) = (alpha*f(int_dofs)  + (1-alpha)*f_prev(int_dofs))* delta_t + mass_mat(int_dofs, int_dofs)*u_0(int_dofs);
    end
    res = f(int_dofs) - alpha * stiff_mat(int_dofs, int_dofs)*u(int_dofs) * delta_t  - mass_mat(int_dofs, int_dofs)*u(int_dofs);
    %jacobian of the residuum
    if (plot_data.print_info)
        fprintf('\n \t Assembly tangent matrix');
    end
    J = alpha * stiff_mat(int_dofs, int_dofs) * delta_t + mass_mat(int_dofs, int_dofs);
    % compute temperature update
    if (plot_data.print_info)
        fprintf('\n \t Evaluate solution increment');
    end
    delta_u = J \ res;
    u(int_dofs) = u(int_dofs) + delta_u;
    if (plot_data.print_info)
        fprintf('\n res norm = %d \n',norm(res));
    end
    
    % STOPPING CRITERIA
    if (iter == problem_data.num_Newton_iter)
        fprintf('\nWarning: Reached max number of iterations = %d \n',iter);
        if (norm(res) < problem_data.Newton_tol)
            problem_data.non_linear_convergence_flag = 1;
        end
        break;
    end
    
end % END NEWTON-RAPHSON ITERATION LOOP

% UPDATE ADAPTIVITY NON-LINEAR CONVERGENCE FLAG
if (norm(res) < problem_data.Newton_tol)
    problem_data.non_linear_convergence_flag = 1;
end

end