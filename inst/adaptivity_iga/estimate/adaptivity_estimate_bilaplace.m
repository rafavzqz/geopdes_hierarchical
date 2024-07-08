% ADAPTIVITY_ESTIMATE_BILAPLACE: Computation of a posteriori error indicators for Bilaplacian problem
%
% USAGE:
%
%   est = adaptivity_estimate_bilaplace (u, hmsh, hspace, problem_data, adaptivity_data)
%
% INPUT:
%
%   u:            degrees of freedom
%   hmsh:         object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace:       object representing the space of hierarchical splines (see hierarchical_space)
%   problem_data: a structure with data of the problem. For this function, it must contain the fields:
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


function est = adaptivity_estimate_bilaplace (u, hmsh, hspace, problem_data, adaptivity_data)

if (isfield(adaptivity_data, 'C0_est'))
    C0_est = adaptivity_data.C0_est;
else
    C0_est = 1;
end

[ders, F] = hspace_eval_hmsh (u, hspace, hmsh, {'bilaplacian'});
der4num = ders{1};

x = cell (hmsh.rdim, 1);
for idim = 1:hmsh.rdim
    x{idim} = reshape (F(idim,:), [], hmsh.nel);
end

if(isfield (problem_data, 'point_load')) %%%%%% be careful ! this works only if I have either a point load or a distributed force
%     pts{1} = problem_data.local_point(1);
%     pts{2} = problem_data.local_point(2);
%     bilaplacian_at_point = sp_eval (u, hspace, problem_data.geometry, pts, 'bilaplacian');
%     
%     rhs = op_point_load_hier(hspace, hmsh, problem_data.point_load, problem_data.local_point);
%     
    global_point = problem_data.geometry.map(problem_data.local_point);
    
    x{1} = x{1} - global_point(1);
    x{2} = x{2} - global_point(2);
    
    d = sqrt(x{1}.^2 + x{2}.^2);
    size_d = size(d);
    rhs = zeros(size_d(1),size_d(2));
    
    for i = 1:size_d(1)
        for j = 1:size_d(2)
            rhs(i,j) = problem_data.approximated_delta_dirac (d(i,j));
        end
    end
%     rhs = problem_data.approximated_delta_dirac (d, 0.08, hmsh.msh_lev{hspace.nlevels}.element_size(1));
    val_D = problem_data.D(x{:});
    aux = (rhs + val_D.*der4num).^2;
else
    valf = problem_data.f (x{:});
    val_D = problem_data.D(x{:});
    aux = (valf + val_D.*der4num).^2; % size(aux) = [hmsh.nqn, hmsh.nel], interior residual at quadrature nodes
end


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
        
        est = sqrt (sum (aux.*w));
        est = C0_est*(h.^2).*est(:);
        
    case 'functions',
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
        disp('Estimate bilaplace for functions is NOT implemented yet !')
end

end