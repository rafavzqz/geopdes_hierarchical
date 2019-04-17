% CHANGE_CONNECTIVITY_LOCALIZED_CSUB: compute the new connectivity related to the localized Csub.
%
% function sp_lev = change_connectivity_localized_Csub (sp_lev, hspace, hmsh)
%
% Compute the new indices and number of dofs related to the localized Csub
%
% INPUT:  
%
%
% OUTPUT:
%
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2017-2019 Rafael Vazquez
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

function sp_lev = change_connectivity_localized_Csub (sp_lev, hspace, hmsh, ilev)

indices = unique (sp_lev.connectivity);
[~,position] = ismember (sp_lev.connectivity, indices);
fun_on_active = sp_get_basis_functions (hspace.space_of_level(ilev), hmsh.mesh_of_level(ilev), hmsh.active{ilev});
fun_on_deact = sp_get_basis_functions (hspace.space_of_level(ilev), hmsh.mesh_of_level(ilev), hmsh.deactivated{ilev});
fun_on_deact = union (fun_on_active, fun_on_deact);
sp_lev.ndof = numel (fun_on_deact);
sp_lev.connectivity = position;

end
