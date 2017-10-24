% ADAPTIVITY_ESTIMATE_POISOON: Computation of a posteriori error indicators for Poisson problem, using globally smooth (C^1) hierarchical spaces.
%
% We consider the diffusion problem
%
%   du(x,t)/dt - div ( epsilon(x) grad (u(x,t))) = f    in Omega = F((0,1)^n)
%                               epsilon(x) du/dn = g    on Gamma_N
%                                              u = h    on Gamma_D
%
% USAGE:
%
%   est = adaptivity_estimate_poison (u, hmsh, hspace, problem_data, adaptivity_data)
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


function est = adaptivity_estimate_poisson (u, u_last, itime, hmsh, hspace, problem_data, adaptivity_data)

if (isfield(adaptivity_data, 'C0_est'))
    C0_est = adaptivity_data.C0_est;
else
    C0_est = 1;
end

%Assume homogeneous time discretization
if itime < numel(problem_data.time_discretization)
    delta_t = problem_data.time_discretization(itime+1) - problem_data.time_discretization(itime);
else
    delta_t = problem_data.time_discretization(itime) - problem_data.time_discretization(itime-1);
end

[ders, F] = hspace_eval_hmsh (u, hspace, hmsh, {'value','gradient', 'laplacian'});
num = ders{1};
dernum = ders{2};
der2num = ders{3};
[ders_last, ~] = hspace_eval_hmsh (u_last, hspace, hmsh, {'value'});
num_last = ders_last{1};

x = cell (hmsh.rdim, 1);
for idim = 1:hmsh.rdim
    x{idim} = reshape (F(idim,:), [], hmsh.nel);
end

aux = 0;
path = cell (hmsh.rdim, 1);
for idim = 1:hmsh.rdim
    path{idim} = ones(size(x{idim})).*problem_data.path(itime, idim);
end
valf = problem_data.f (x{:}, path{:});
val_c_diff = problem_data.c_diff(x{:});
val_c_cap = problem_data.c_cap(x{:});
if (isfield (problem_data, 'grad_c_diff'))
    val_grad_c_diff  = feval (problem_data.grad_c_diff, x{:});
    aux = reshape (sum (val_grad_c_diff .* dernum, 1), size(valf));
end
aux = (valf + val_c_diff.*der2num + aux - val_c_cap.*(num - num_last)./delta_t).^2; % size(aux) = [hmsh.nqn, hmsh.nel], interior residual at quadrature nodes

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

end