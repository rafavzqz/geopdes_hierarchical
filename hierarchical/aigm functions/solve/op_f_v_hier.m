% OP_F_V_HIER: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i) for hierarchical splines, 
%  exploiting the multilevel structure.
%
%   rhs = op_f_v_hier (hspace, hmsh, coeff);
%
% INPUT:
%     
%   hspace: object representing the hierarchical space of trial functions (see XXXXXXXXX)
%   hmsh:   object representing the hierarchical mesh (see XXXXXXXXXXXXX)
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

function rhs = op_f_v_hier (hspace, hmsh, f)

  rhs = zeros (hspace.ndof, 1);

  ndof_per_level = hspace.ndof_per_level;
  dif = hmsh.nlevels - hspace.nlevels;
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%  CAN WE REMOVE THE WARNING (and dif)?
  if (dif)
    warning('Different number of levels for mesh and space. This should never happen')
    ndof_per_level = [ndof_per_level(:); zeros(dif,1)];
  end
  
  ndofs = 0;
  for ilev = 1:hmsh.nlevels
    ndofs = ndofs + ndof_per_level(ilev); %hspace.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)
      x = cell (hmsh.rdim, 1);
      for idim = 1:hmsh.rdim
        x{idim} = hmsh.msh_lev{ilev}.geo_map(idim,:,:);
      end
      b_lev = op_f_v (hspace.sp_lev{ilev}, hmsh.msh_lev{ilev}, f(x{:}));

      dofs = 1:ndofs;
      rhs(dofs) = rhs(dofs) + hspace.C{ilev}'*b_lev;
    end
  end

end
