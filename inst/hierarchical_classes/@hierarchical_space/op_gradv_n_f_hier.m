% OP_F_V_HIER: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, (grad v n)_i) for hierarchical splines, 
%  exploiting the multilevel structure.
%
%   rhs = op_gradv_n_f_hier (hspace, hmsh, coeff);
%
% INPUT:
%     
%   hspace: object representing the hierarchical space of test functions (see hierarchical_space)
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   coeff:  function handle to compute the source function
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

function rhs = op_gradv_n_f_hier (hspace, hmsh, f, iside)

rhs = zeros (hspace.ndof, 1);

ndofs = 0;

hmsh_sfi = hmsh_boundary_side_from_interior (hmsh, iside);

x = cell (hmsh.rdim, 1);
shifting_indices = cumsum ([0 hmsh.boundary(iside).nel_per_level]);
for ilev = 1:hmsh.boundary(iside).nlevels
    ndofs = ndofs + hspace.ndof_per_level(ilev);

    if (hmsh.boundary(iside).nel_per_level(ilev) > 0)
      elements = shifting_indices(ilev)+1:shifting_indices(ilev+1);
      msh_side = hmsh_eval_boundary_side (hmsh, iside, elements);
      msh_side_from_interior = hmsh_sfi.mesh_of_level(ilev);

      sp_bnd = hspace.space_of_level(ilev).constructor (msh_side_from_interior);
      msh_side_from_interior_struct = msh_evaluate_element_list (msh_side_from_interior, hmsh_sfi.active{ilev});
      sp_bnd_struct = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', true, 'gradient', true);

      for idim = 1:hmsh.rdim
        x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
      end
      sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true); %%%% find a smarter way here
      indices = unique (sp_lev.connectivity);
      [~,position] = ismember (sp_bnd_struct.connectivity, indices);
      fun_on_active = sp_get_basis_functions (hspace.space_of_level(ilev), hmsh.mesh_of_level(ilev), hmsh.active{ilev});
      fun_on_deact = sp_get_basis_functions (hspace.space_of_level(ilev), hmsh.mesh_of_level(ilev), hmsh.deactivated{ilev});
      fun_on_deact = union (fun_on_active, fun_on_deact);
      sp_bnd_struct.ndof = numel (fun_on_deact);
      sp_bnd_struct.connectivity = position;

      b_lev = op_gradv_n_f (sp_bnd_struct, msh_side, f(x{:}, iside));

      dofs = 1:ndofs;
      rhs(dofs) = rhs(dofs) + hspace.Csub{ilev}.' * b_lev;
    end
end
  
%   rhs = zeros (hspace.ndof, 1);
% 
%   ndofs = 0;
%   for ilev = 1:hmsh.nlevels
%     ndofs = ndofs + hspace.ndof_per_level(ilev);
%     if (ilev <= numel(hmsh.boundary(iside).nel_per_level) && hmsh.boundary(iside).nel_per_level(ilev) > 0)
%       x = cell (hmsh.rdim, 1);
%       %%% use this trick to store the boundary information on the fly
%       %%% inside mesh_of_level
%       hmsh.mesh_of_level(ilev).boundary = hmsh.boundary(iside).mesh_of_level(ilev);
%       msh_side = msh_eval_boundary_side (hmsh.mesh_of_level(ilev), iside, hmsh.boundary(iside).active{ilev});
%       msh_side_from_interior = msh_boundary_side_from_interior (hmsh.mesh_of_level(ilev), iside);
% 
%       sp_bnd = hspace.space_of_level(ilev).constructor (msh_side_from_interior);
%       msh_side_from_interior_struct = msh_evaluate_element_list(msh_side_from_interior, hmsh.boundary(iside).active{ilev});
%       sp_bnd = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', true, 'gradient', true);
% 
%       for idim = 1:hmsh.rdim
%         x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
%       end
% 
%       b_lev = op_gradv_n_f (sp_bnd, msh_side, f(x{:}, iside));
% 
%       dofs = 1:ndofs;
%       rhs(dofs) = rhs(dofs) + hspace.Csub{ilev}.' * b_lev;
%     end
%   end

end
