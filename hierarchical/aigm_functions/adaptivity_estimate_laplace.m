% ADAPTIVITY_ESTIMATE_LAPLACE: Computation of a posteriori error indicators for Laplacian problem, using smooth (C^1) hierarchical spaces.
%
% We consider the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^n)
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% USAGE:
%
% est = adaptivity_estimate_laplace (u, hmsh, hspace, problem_data, adaptivity_data)
%
% INPUT:
%
%   u:      degrees of freedom
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the space of hierarchical splines (see hierarchical_space)
%   problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - c_diff:       diffusion coefficient (epsilon in the equation)
%    - grad_c_diff:  gradient of the diffusion coefficient 
%    - f:            source term
%   adaptivity_data: a structure with the data for the adaptivity method. In particular, it contains the fields:
%    - flag:         'elements' or 'functions', depending on the refinement strategy.
%    - C0_est      : multiplicative constant for the error indicators 
%                    
%
% OUTPUT:
%
%   est: computed a posteriori error indicators
%           - (Buffa and Giannelli, 2016) When adaptivity_data.flag == 'elements': for an element Q, 
%                          est_Q := C0_est*h_Q*(int_Q |f + div(epsilon(x) grad(U))|^2)^(1/2), 
%           where h_Q is the local meshsize and U is the Galerkin solution
%           - (Buffa and Garau, 2016) When adaptivity_data.flag == 'functions': for a B-spline basis function b, 
%                          est_b := C0_est*h_b*(int_{supp b} a_b*|f + div(epsilon(x) grad(U))|^2*b)^(1/2), 
%           where h_b is the local meshsize, a_b is the coefficient of b for the partition-of-unity in the hierarchical basis, and U is the Galerkin solution
%
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

% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% - When c_diff is constant, we could use a flag in order to avoid some useless computations
% - Up to now, we are not considering jumps; thus, we assume also that c_diff is smooth. 
% - These a posteriori error indicators a designed for homogeneous boundary conditions. I will think how to modify them for non-homogeneous boundary conditions
%

function est = adaptivity_estimate_laplace (u, hmsh, hspace, problem_data, adaptivity_data)

[der2num, F] = hspace_eval_hmsh (u, hspace, hmsh, 'laplacian'); 
[dernum, ~] = hspace_eval_hmsh (u, hspace, hmsh, 'gradient'); 


x = cell (hmsh.rdim, 1);
for idim = 1:hmsh.rdim;
  x{idim} = reshape (F(idim,:), [], hmsh.nel);
end

val_c_diff = problem_data.c_diff(x{:});
val_grad_c_diff  = feval (problem_data.grad_c_diff, x{:});

valf = problem_data.f (x{:});

aux = (valf + val_c_diff.*der2num + squeeze(dot(val_grad_c_diff,dernum,1))).^2; % size(aux) = [hmsh.nqn, hmsh.nel], interior residual at quadrature nodes

quad_weights = [];
jacdet = [];
h = [];
ms = zeros (hmsh.nlevels, 1);
for ilev = 1:hmsh.nlevels % Active levels
    if (hmsh.msh_lev{ilev}.nel ~= 0)
            quad_weights = cat(2,quad_weights, hmsh.msh_lev{ilev}.quad_weights);
            jacdet = cat(2,jacdet, hmsh.msh_lev{ilev}.jacdet);
            h = cat (1, h, hmsh.msh_lev{ilev}.element_size(:));
            ms(ilev) = max (hmsh.msh_lev{ilev}.element_size);
    else
        ms(ilev) = 0;
    end
end
w = quad_weights .* jacdet;
h = h * sqrt (hmsh.ndim);
ms = ms * sqrt (hmsh.ndim);

switch adaptivity_data.flag
    case 'elements',
        est = sqrt (sum (aux.*w));
        est = adaptivity_data.C0_est*h.*est(:);
        
    case 'functions',
        Nf = cumsum ([0; hspace.ndof_per_level(:)]);
        dof_level = zeros (hspace.ndof, 1);
        for lev = 1:hspace.nlevels
            dof_level(Nf(lev)+1:Nf(lev+1)) = lev;
        end
        coef = ms(dof_level).*sqrt(hspace.coeff_pou(:));
        
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        % If we allow that the third input of op_f_v_hier can be f at the
        % quadrature nodes, we can use the following line instead of the
        % lines below.
        % est = adaptivity_data.C0_est*coef.* sqrt(op_f_v_hier (hspace, hmsh, aux));
        
        est = zeros(hspace.ndof,1);
        ndofs = 0;
        Ne = cumsum([0; hmsh.nel_per_level(:)]);
        for ilev = 1:hmsh.nlevels 
            ndofs = ndofs + hspace.ndof_per_level(ilev);
            if (hmsh.nel_per_level(ilev) > 0)
                ind_e = (Ne(ilev)+1):Ne(ilev+1);
                sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
                b_lev = op_f_v (sp_lev, hmsh.msh_lev{ilev}, aux(:,ind_e));
                dofs = 1:ndofs; 
                est(dofs) = est(dofs) + hspace.C{ilev}'*b_lev;
            end
        end
        est = adaptivity_data.C0_est*coef.*sqrt(est);
end
