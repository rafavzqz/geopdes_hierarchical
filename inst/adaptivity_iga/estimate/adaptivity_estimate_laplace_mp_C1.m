% ADAPTIVITY_ESTIMATE_LAPLACE: Computation of a posteriori error indicators for Laplacian problem, using globally smooth (C^1) hierarchical spaces.
%
% We consider the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^n)
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% USAGE:
%
%   est = adaptivity_estimate_laplace (u, hmsh, hspace, problem_data, adaptivity_data)
%
% INPUT:
%
%   u:            degrees of freedom
%   hmsh:         object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace:       object representing the space of hierarchical splines (see hierarchical_space)
%   problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - c_diff:        diffusion coefficient (epsilon in the equation), assumed to be a smooth (C^1) function
%    - grad_c_diff:   gradient of the diffusion coefficient (equal to zero if not present)
%    - f:             function handle of the source term
%   adaptivity_data: a structure with the data for the adaptivity method. In particular, it contains the fields:
%    - flag:          'elements' or 'functions', depending on the refinement strategy.
%    - C0_est:        multiplicative constant for the error indicators
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
%           where h_b is the local meshsize, a_b is the coefficient of b for the partition-of-unity, and U is the Galerkin solution
%
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
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


function est = adaptivity_estimate_laplace_mp_C1 (u, hmsh, hspace, problem_data, adaptivity_data)

if (isfield(adaptivity_data, 'C0_est'))
    C0_est = adaptivity_data.C0_est;
else
    C0_est = 1;
end

[ders, F] = hspace_eval_hmsh (u, hspace, hmsh, {'gradient', 'laplacian'});
dernum = ders{1};
der2num = ders{2};

x = cell (hmsh.rdim, 1);
for idim = 1:hmsh.rdim
    x{idim} = reshape (F(idim,:), [], hmsh.nel);
end

aux = 0;
valf = problem_data.f (x{:});
val_c_diff = problem_data.c_diff(x{:});
if (isfield (problem_data, 'grad_c_diff'))
    val_grad_c_diff  = feval (problem_data.grad_c_diff, x{:});
    aux = reshape (sum (val_grad_c_diff .* dernum, 1), size(valf));
end
aux = (valf + val_c_diff.*der2num + aux).^2; % size(aux) = [hmsh.nqn, hmsh.nel], interior residual at quadrature nodes

switch adaptivity_data.flag
    case 'elements'
        w = [];
        h = [];
        for ilev = 1:hmsh.nlevels
            if (hmsh.msh_lev{ilev}.nel ~= 0)
                w = cat (2, w, hmsh.msh_lev{ilev}.quad_weights .* hmsh.msh_lev{ilev}.jacdet);
                h = cat (1, h, hmsh.msh_lev{ilev}.element_size(:));
            end
        end
        h = h * sqrt (hmsh.ndim);
        
        est = sqrt (sum (aux.*w));
        est = C0_est*h.*est(:);
        
    case 'functions'
        error ('The estimator by functions is not implemented for the C^1 case')
%         ms = zeros (hmsh.nlevels, 1);
%         for ilev = 1:hmsh.nlevels
%             if (hmsh.msh_lev{ilev}.nel ~= 0)
%                 ms(ilev) = max (hmsh.msh_lev{ilev}.element_size);
%             else
%                 ms(ilev) = 0;
%             end
%         end
%         ms = ms * sqrt (hmsh.ndim);
%         
%         Nf = cumsum ([0; hspace.ndof_per_level(:)]);
%         dof_level = zeros (hspace.ndof, 1);
%         for lev = 1:hspace.nlevels
%             dof_level(Nf(lev)+1:Nf(lev+1)) = lev;
%         end
%         coef = ms(dof_level).*sqrt(hspace.coeff_pou(:));
%         
%         est = zeros(hspace.ndof,1);
%         ndofs = 0;
%         Ne = cumsum([0; hmsh.nel_per_level(:)]);
%         for ilev = 1:hmsh.nlevels
%             ndofs = ndofs + hspace.ndof_per_level(ilev);
%             if (hmsh.nel_per_level(ilev) > 0)
%                 ind_e = (Ne(ilev)+1):Ne(ilev+1);
%                 sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
%                 b_lev = op_f_v (sp_lev, hmsh.msh_lev{ilev}, aux(:,ind_e));
%                 dofs = 1:ndofs;
%                 est(dofs) = est(dofs) + hspace.Csub{ilev}.' * b_lev;
%             end
%         end
%         est = C0_est * coef .* sqrt(est);
end

end