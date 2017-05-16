% OP_F_V_TIME_HIER: assemble the right-hand side time-dependent vector r(t) = [r(i)], with  r(i) = (f(t), v_i) for hierarchical splines,
%  exploiting the multilevel structure.
%
%   rhs = op_f_v_time_hier (hspace, hmsh, problem_data);
%
% INPUT:
%
%   hspace: object representing the hierarchical space of test functions (see hierarchical_space_mp)
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh_mp)
%   problem_data: a structure with data of the problem. It contains the fields:
%    - f:                   function handle of the source term
%    - path:                position of the source term
%    - time_discretization: time value at each time step
%   time_step:              actual time step
%
% OUTPUT:
%
%   rhs: assembled right-hand side
%
% The multilevel structure is exploited in such a way that only functions
%  of the same level of the active elements have to be computed. See also
%  op_gradu_gradv_hier for more details.
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

function rhs = op_f_v_time_hier (hspace, hmsh, problem_data, time_step)

rhs = zeros (hspace.ndof, 1);
ndofs = 0;

for ilev = 1:hmsh.nlevels
    ndofs = ndofs + hspace.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)
        x = cell (hmsh.rdim, 1);
        path = cell (hmsh.rdim, 1);
        for idim = 1:hmsh.rdim
            x{idim} = reshape (hmsh.msh_lev{ilev}.geo_map(idim,:,:), hmsh.mesh_of_level(ilev).nqn, hmsh.nel_per_level(ilev));
            path{idim} = ones(size(x{idim})).*problem_data.path(time_step, idim);
        end
        sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
        b_lev = op_f_v (sp_lev, hmsh.msh_lev{ilev}, problem_data.f(x{:}, path{:}));
        
        dofs = 1:ndofs;
        rhs(dofs) = rhs(dofs) + hspace.Csub{ilev}.' * b_lev;
    end
end

end
