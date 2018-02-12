% ADAPTIVITY_ESTIMATE_LINEAR_EL_mixed: Computation of a posteriori error indicators for linear elasticity problem, using globally smooth (C^1) hierarchical spaces,
% in the case of mixed formulation.
%
% NOTE: works for 'elements' case
%
% We consider the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^n)
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% USAGE:
%
%   est = adaptivity_estimate_lINEAR_EL (u, hmsh, hspace, problem_data, adaptivity_data)
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
% Copyright (C) 2017 Cesare Bracco, Rafael Vazquez
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


function est = adaptivity_estimate_linear_el_mixed (u, press, hmsh, hspace_u, hspace_p, problem_data, adaptivity_data)

if (isfield(adaptivity_data, 'C0_est'))
    C0_est = adaptivity_data.C0_est;
else
    C0_est = 1;
end

[ders2_u, F] = hspace_eval_hmsh (u, hspace_u, hmsh, 'hessian');
[divergence_u, F] = hspace_eval_hmsh (u, hspace_u, hmsh, 'divergence');
[gradient_p, F] = hspace_eval_hmsh (press, hspace_p, hmsh, 'gradient');
[pressure, F] = hspace_eval_hmsh (press, hspace_p, hmsh);

x = cell (hmsh.rdim, 1);
for idim = 1:hmsh.rdim;  %rdim is the dimension of the physical domain
    x{idim} = reshape (F(idim,:), [], hmsh.nel);
end

lambda=problem_data.lambda_lame(1); %we assume that lambda and mu are constants
mu=problem_data.mu_lame(1);

valf = problem_data.f (x{1},x{2});
for h = 1:hspace_u.ncomp
    partials_a = 0;
    partials_b = 0;
    for t=1:hspace_u.ncomp
        if t~=h
            partials_a = partials_a + squeeze(ders2_u(t,h,t,:,:)); %mixed derivatives of all components (except h-th component)
            partials_b = partials_b + squeeze(ders2_u(h,t,t,:,:)); %second derivatives of h-th component (except w.r.t. h-th variable)
        end
    end
    
%     size(gradient_p(h,:,:))
%     size((2*mu)*squeeze(ders2_u(h,h,h,:,:)) + (mu)*partials_a + mu*partials_b)
    divergence(h,:,:) = (2*mu)*squeeze(ders2_u(h,h,h,:,:)) + (mu)*partials_a + mu*partials_b + squeeze(gradient_p(h,:,:));
end
aux1=(valf+divergence).^2;  %residual 1
aux1=squeeze(sum(aux1));

% size(press)
% size(divergence_u)

aux2=((1/lambda)*pressure-divergence_u).^2;   %residual 2
%aux2=squeeze(sum(aux2));

switch adaptivity_data.flag
    case 'elements',
        w = [];
        h = [];
        for ilev = 1:hmsh.nlevels
            if (hmsh.msh_lev{ilev}.nel ~= 0)
                w = cat (2, w, hmsh.msh_lev{ilev}.quad_weights .* hmsh.msh_lev{ilev}.jacdet);
                h = cat (1, h, hmsh.msh_lev{ilev}.element_size(:));
            end
        end
        h = h * sqrt (hmsh.ndim);
        
%         size(aux1)
%         size(aux2)
        
        est1 = sum (aux1.*w);
        est2 = sum (aux2.*w);
        est = C0_est*sqrt((h.^2).*est1(:)+est2(:));
        
    case 'functions',
        ms = zeros (hmsh.nlevels, 1);
        for ilev = 1:hmsh.nlevels
            if (hmsh.msh_lev{ilev}.nel ~= 0)
                ms(ilev) = max (hmsh.msh_lev{ilev}.element_size);
            else
                ms(ilev) = 0;
            end
        end
        ms = ms * sqrt (hmsh.ndim);
        
        Nf = cumsum ([0; hspace.ndof_per_level(:)]);
        dof_level = zeros (hspace.ndof, 1);
        for lev = 1:hspace.nlevels
            dof_level(Nf(lev)+1:Nf(lev+1)) = lev;
        end
        coef = ms(dof_level).*sqrt(hspace.coeff_pou(:));
        
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
                est(dofs) = est(dofs) + hspace.Csub{ilev}.' * b_lev;
            end
        end
        est = C0_est * coef .* sqrt(est);
end

end