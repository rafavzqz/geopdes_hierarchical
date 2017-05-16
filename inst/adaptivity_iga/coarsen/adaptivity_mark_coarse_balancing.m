%ADAPTIVITY_MARK_COARSE_BALANCING relax the mesh avoiding two level jumps at the element
%interface
%
%     [marked_elements] = adaptivity_mark_coarse_balancing (hmsh, marked_elements, ringIndex)
%
%
% INPUT:
%
%     hmsh: the hierarchical mesh (see hierarchical_mesh)
%     marked_elements: marked elements from previous marking strategy
%     ringIndex: max neighbours to be check in balancing
%
% OUTPUT:
%
%     marked_elements: marked elements to relax the hierarchical mesh
%
% Copyright (C) 2017 Massimo Carraturo
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
function  [marked_elements]  = adaptivity_mark_coarse_balancing(hmsh, marked_elements, ringIndex)

active = hmsh.active;
el_dir = zeros(3, 1);
%loop over level
for lev=hmsh.nlevels:-1:2
    %loop over active elements of the level lev
    for el=1:numel(active{lev})
        %get parents of the element
        [parent, ~] = hmsh_get_parent(hmsh, lev, active{lev}(el));
        %get neighbours of the parent
        % first get the local indeces in each direction ...
            [el_dir(1), el_dir(2), el_dir(3)] = ind2sub(hmsh.mesh_of_level(lev-1).nel_dir, parent);
        % ... then get the parent's neighbours
            neighbours  = get_neighbours(hmsh, lev-1, el_dir(1), el_dir(2), el_dir(3), ringIndex);
        %get active neighbours
            active_neighbours = intersect( neighbours, active{lev-1});
        %check if active neighbours are on the marked element list
        if (~isempty(intersect(active_neighbours, marked_elements{lev-1})))
            %remove them from the list
            marked_elements{lev-1} = setdiff(marked_elements{lev-1},active_neighbours);
        end
    end
    
end

end

function  neighbours  = get_neighbours( hmsh, lev, ind_x, ind_y, ind_z, relaxataion_par )
% x direction
neighbours_x_dir = linspace(ind_x-relaxataion_par, ind_x+relaxataion_par,2*relaxataion_par+1);
neighbours_x = zeros(size(neighbours_x_dir));
for i=1:numel(neighbours_x_dir)
    if (neighbours_x_dir(i)<=hmsh.mesh_of_level(lev).nel_dir(1) && neighbours_x_dir(i)>0)
        neighbours_x(i) = sub2ind(hmsh.mesh_of_level(lev).nel_dir, neighbours_x_dir(i), ind_y, ind_z);
    end
end
neighbours_x = unique(neighbours_x(neighbours_x>0));
% y direction
neighbours_y = [];
if numel(hmsh.mesh_of_level(lev).nel_dir) > 1
    neighbours_y_dir = linspace(ind_y-relaxataion_par, ind_y+relaxataion_par,2*relaxataion_par+1);
    neighbours_y = zeros(size(neighbours_y_dir));
    for i=1:numel(neighbours_y_dir)
        if(neighbours_y_dir(i)<=hmsh.mesh_of_level(lev).nel_dir(2)  && neighbours_y_dir(i)>0)
            neighbours_y(i) = sub2ind(hmsh.mesh_of_level(lev).nel_dir, ind_x, neighbours_y_dir(i), ind_z);
        end
    end
end
neighbours_y = unique(neighbours_y(neighbours_y>0));
% z direction
neighbours_z = [];
if numel(hmsh.mesh_of_level(lev).nel_dir) > 2
    neighbours_z_dir = linspace(ind_z-relaxataion_par, ind_z+relaxataion_par,2*relaxataion_par+1);
    neighbours_z = zeros(size(neighbours_z_dir));
    for i=1:numel(neighbours_z_dir)
        if (neighbours_z_dir(i)<=hmsh.mesh_of_level(lev).nel_dir(3)  && neighbours_z_dir(i)>0)
            neighbours_z(i) = sub2ind(hmsh.mesh_of_level(lev).nel_dir, ind_x, ind_y, neighbours_z_dir(i));% fill an auxiliary list of elements
        end
    end
end
neighbours_z = unique(neighbours_z(neighbours_z>0));

neighbours = unique([neighbours_x neighbours_y neighbours_z]);
end